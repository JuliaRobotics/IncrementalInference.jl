# clique state machine for tree based initialization and inference

# newer exports
# export towardUpOrDwnSolve_StateMachine, checkIfCliqNullBlock_StateMachine, doAnyChildrenNeedDwn_StateMachine
# export mustInitUpCliq_StateMachine, doCliqUpSolveInitialized_StateMachine
# export rmUpLikeliSaveSubFg_StateMachine
# export blockCliqSiblingsParentChildrenNeedDown_StateMachine

export  doCliqDownSolve_StateMachine,
        cleanupAfterDownSolve_StateMachine,
        specialCaseRootDownSolve_StateMachine,
        canCliqDownSolve_StateMachine,
        finishCliqSolveCheck_StateMachine,
        mustInitUpCliq_StateMachine,
        doCliqUpSolveInitialized_StateMachine,
        rmUpLikeliSaveSubFg_StateMachine,
        somebodyLovesMe_StateMachine,
        waitChangeOnParentCondition_StateMachine,
        slowOnPrntAsChildrNeedDwn_StateMachine,
        towardUpOrDwnSolve_StateMachine,
        canCliqMargSkipUpSolve_StateMachine,
        attemptDownInit_StateMachine,
        rmMsgLikelihoodsAfterDwn_StateMachine,
        blockUntilSiblingsStatus_StateMachine,
        slowIfChildrenNotUpSolved_StateMachine,
        blockUntilChildrenHaveStatus_StateMachine,
        dwnInitSiblingWaitOrder_StateMachine,
        collectDwnInitMsgFromParent_StateMachine,
        trafficRedirectConsolidate459_StateMachine,
        doAllSiblingsNeedDwn_StateMachine,
        checkIfCliqNullBlock_StateMachine,
        determineCliqNeedDownMsg_StateMachine,
        doAnyChildrenNeedDwn_StateMachine,
        decideUpMsgOrInit_StateMachine,
        doesParentNeedDwn_StateMachine,
        attemptCliqInitUp_StateMachine,
        sendCurrentUpMsg_StateMachine,
        buildCliqSubgraph_StateMachine,
        buildCliqSubgraphForDown_StateMachine,
        isCliqUpSolved_StateMachine,
        checkChildrenAllUpRecycled_StateMachine,
        testCliqCanIncremtUpdate_StateMachine!,
        testCliqCanRecycled_StateMachine



"""
    $SIGNATURES

Specialized info logger print function to show clique state machine information
in a standardized form.
"""
function infocsm(csmc::CliqStateMachineContainer, str::A) where {A <: AbstractString}

  tm = string(Dates.now())
  tmt = split(tm, 'T')[end]

  lbl = getLabel(csmc.cliq)
  lbl1 = split(lbl,',')[1]
  cliqst = getCliqueStatus(csmc.cliq)

  with_logger(csmc.logger) do
    @info "$tmt | $(csmc.cliq.index)---$lbl1 @ $(cliqst) | "*str
  end
  flush(csmc.logger.stream)
  nothing
end

"""
    $SIGNATURES

Do cliq downward inference

Notes:
- State machine function nr. 11
"""
function doCliqDownSolve_StateMachine(csmc::CliqStateMachineContainer)
  infocsm(csmc, "11, doCliqDownSolve_StateMachine")
  setCliqDrawColor(csmc.cliq, "red")
  opts = getSolverParams(csmc.dfg)

  # get down msg from parent (assuming root clique CSM wont make it here)
  # this looks like a pull model #674
  prnt = getParent(csmc.tree, csmc.cliq)
  dwnmsgs = fetchMsgDwnThis(prnt[1])
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- dwnmsgs=$(collect(keys(dwnmsgs.belief)))")

  # maybe cycle through separators (or better yet, just use values directly -- see next line)
  msgfcts = addMsgFactors!(csmc.cliqSubFg, dwnmsgs, DownwardPass)
  # force separator variables in cliqSubFg to adopt down message values
  updateSubFgFromDownMsgs!(csmc.cliqSubFg, dwnmsgs, getCliqSeparatorVarIds(csmc.cliq))

  # add required all frontal connected factors
  if !opts.useMsgLikelihoods
    newvars, newfcts = addDownVariableFactors!(csmc.dfg, csmc.cliqSubFg, csmc.cliq, csmc.logger, solvable=1)
  end

  # store the cliqSubFg for later debugging
  if opts.dbg
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforedownsolve"))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforedownsolve.pdf"))
  end

  ## new way
  # calculate belief on each of the frontal variables and iterate if required
  solveCliqDownFrontalProducts!(csmc.cliqSubFg, csmc.cliq, opts, csmc.logger)

  # compute new down messages
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- going to set new down msgs.")
  newDwnMsgs = getSetDownMessagesComplete!(csmc.cliqSubFg, csmc.cliq, dwnmsgs, csmc.logger)
  
    #JT 459 putMsgDwnThis!(cliq, newDwnMsgs), DF still looks like a pull model here #674
    putMsgDwnThis!(csmc.cliq.data, newDwnMsgs, from=:putMsgDwnThis!) # putCliqueMsgDown!

    # update clique subgraph with new status
    setCliqDrawColor(csmc.cliq, "lightblue")

  csmc.dodownsolve = false
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- finished with downGibbsCliqueDensity, now update csmc")

  # go to 11b.
  return cleanupAfterDownSolve_StateMachine
end

"""
    $SIGNATURES

One of the last steps in CSM to clean up after a down solve.

Notes
- CSM function 11b.
"""
function cleanupAfterDownSolve_StateMachine(csmc::CliqStateMachineContainer)
  # RECENT split from 11 (using #760 solution for deleteMsgFactors)
  opts = getSolverParams(csmc.cliqSubFg)

  # set PPE and solved for all frontals
  for sym in getCliqFrontalVarIds(csmc.cliq)
    # set PPE in cliqSubFg
    setVariablePosteriorEstimates!(csmc.cliqSubFg, sym)
    # set solved flag
    vari = getVariable(csmc.cliqSubFg, sym)
    setSolvedCount!(vari, getSolvedCount(vari, :default)+1, :default )
  end

  # store the cliqSubFg for later debugging
  if opts.dbg
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_afterdownsolve"))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_afterdownsolve.pdf"))
  end

  # transfer results to main factor graph
  frsyms = getCliqFrontalVarIds(csmc.cliq)
  infocsm(csmc, "11, finishingCliq -- going for transferUpdateSubGraph! on $frsyms")
  transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frsyms, csmc.logger, updatePPE=true)

  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- before notifyCliqDownInitStatus!")
  notifyCliqDownInitStatus!(csmc.cliq, :downsolved, logger=csmc.logger)
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- just notified notifyCliqDownInitStatus!")

  # remove msg factors that were added to the subfg
  rmFcts = lsf(csmc.cliqSubFg, tags=[:LIKELIHOODMESSAGE;]) .|> x -> getFactor(csmc.cliqSubFg, x)
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- removing all up/dwn message factors, length=$(length(rmFcts))")
  deleteMsgFactors!(csmc.cliqSubFg, rmFcts) # msgfcts # TODO, use tags=[:LIKELIHOODMESSAGE], see #760

  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- finished, exiting CSM on clique=$(csmc.cliq.index)")
  # and finished
  return IncrementalInference.exitStateMachine
end


"""
$SIGNATURES

Root clique upsolve and downsolve are equivalent, so skip a repeat downsolve, just set messages and just exit directly.

Notes
- State machine function nr. 10b
- Separate out during #459 dwnMsg consolidation.

DevNotes
- TODO should this consolidate some work with 11b?
"""
function specialCaseRootDownSolve_StateMachine(csmc::CliqStateMachineContainer)
  # this is the root clique, so assume already downsolved -- only special case
  dwnmsgs = getCliqDownMsgsAfterDownSolve(csmc.cliqSubFg, csmc.cliq)
  setCliqDrawColor(csmc.cliq, "lightblue")

  # this part looks like a pull model
  # JT 459 putMsgDwnThis!(csmc.cliq, dwnmsgs)
  putMsgDwnThis!(csmc.cliq.data, dwnmsgs, from=:putMsgDwnThis!) # putCliqueMsgDown!
  setCliqueStatus!(csmc.cliq, :downsolved)
  csmc.dodownsolve = false

  # Update estimates and transfer back to the graph
  frsyms = getCliqFrontalVarIds(csmc.cliq)
  # set PPE and solved for all frontals
  for sym in frsyms
    # set PPE in cliqSubFg
    setVariablePosteriorEstimates!(csmc.cliqSubFg, sym)
    # set solved flag
    vari = getVariable(csmc.cliqSubFg, sym)
    setSolvedCount!(vari, getSolvedCount(vari, :default)+1, :default )
  end
  # Transfer to parent graph
  transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frsyms, updatePPE=true)

  notifyCliqDownInitStatus!(csmc.cliq, :downsolved, logger=csmc.logger)

  # bye
  return IncrementalInference.exitStateMachine
end



"""
    $SIGNATURES

Direct state machine to continue with downward solve or exit.

Notes
- State machine function nr. 10
"""
function canCliqDownSolve_StateMachine(csmc::CliqStateMachineContainer)
  infocsm(csmc, "10, canCliqDownSolve_StateMachine, csmc.dodownsolve=$(csmc.dodownsolve).")
  # finished and exit downsolve
  if !csmc.dodownsolve
    infocsm(csmc, "10, canCliqDownSolve_StateMachine -- shortcut exit since downsolve not required.")
    return IncrementalInference.exitStateMachine
  end

  # assume separate down solve via solveCliq! call, but need a csmc.cliqSubFg
  # could be dedicated downsolve that was skipped during previous upsolve only call
  # e.g. federated solving case (or debug)
  if length(ls(csmc.cliqSubFg)) == 0
    # first need to fetch cliq sub graph
    infocsm(csmc, "10, canCliqDownSolve_StateMachine, oops no cliqSubFg detected, lets go fetch a copy first.")

    # go to 2b
    return buildCliqSubgraphForDown_StateMachine
  end
  
  # both parent or otherwise might start by immediately doing downsolve, so likely need cliqSubFg in both cases
  # e.g. federated solving case (or debug)
  prnt = getParent(csmc.tree, csmc.cliq)
  if 0 == length(prnt) # check if have parent
    # go to 10b
    return specialCaseRootDownSolve_StateMachine
  end
  
  # go to 10a
  return somebodyLovesMe_StateMachine
  
end


"""
    $SIGNATURES

Need description for this???

Notes
- State machine function nr.9
"""
function finishCliqSolveCheck_StateMachine(csmc::CliqStateMachineContainer)
  cliqst = getCliqueStatus(csmc.cliq)
  infocsm(csmc, "9, finishCliqSolveCheck_StateMachine")
  if cliqst == :upsolved
      frsyms = getCliqFrontalVarIds(csmc.cliq)
    infocsm(csmc, "9, finishCliqSolveCheck_StateMachine -- going for transferUpdateSubGraph! on $frsyms")
    # TODO what about down solve??
    transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frsyms, csmc.logger, updatePPE=false)

    # remove any solvable upward cached data -- TODO will have to be changed for long down partial chains
    # assuming maximally complte up solved cliq at this point
    lockUpStatus!(csmc.cliq, csmc.cliq.index, true, csmc.logger, true, "9.finishCliqSolveCheck")
    sdims = Dict{Symbol,Float64}()
    for varid in getCliqAllVarIds(csmc.cliq)
      sdims[varid] = 0.0
    end
    updateCliqSolvableDims!(csmc.cliq, sdims, csmc.logger)
    unlockUpStatus!(csmc.cliq)

    # go to 10
    return canCliqDownSolve_StateMachine # IncrementalInference.exitStateMachine
  elseif cliqst == :initialized
    setCliqDrawColor(csmc.cliq, "sienna")

    # go to 7
    return determineCliqNeedDownMsg_StateMachine
  else
    infocsm(csmc, "9, finishCliqSolveCheck_StateMachine -- init not complete and should wait on init down message.")
    setCliqDrawColor(csmc.cliq, "coral")
    # TODO, potential problem with trying to downsolve
    # return canCliqMargSkipUpSolve_StateMachine
  end

  # go to 4b (redirected here during #459 dwnMsg effort)
  return trafficRedirectConsolidate459_StateMachine
    # # go to 4
    # return canCliqMargSkipUpSolve_StateMachine # whileCliqNotSolved_StateMachine
end


"""
    $SIGNATURES

Do up initialization calculations, loosely translates to solving Chapman-Kolmogorov
transit integral in upward direction.

Notes
- State machine function nr. 8f
- Includes initialization routines.
- Adds LIKELIHOODMESSAGE factors but does not remove.
- gets msg likelihoods from cliqSubFg, see #760

DevNotes
- TODO: Make multi-core
"""
function mustInitUpCliq_StateMachine(csmc::CliqStateMachineContainer)
  setCliqDrawColor(csmc.cliq, "green")

  # check if init is required and possible
  infocsm(csmc, "8f, mustInitUpCliq_StateMachine -- going for doCliqAutoInitUpPart1!.")
  # get incoming clique up messages
  upmsgs = getMsgsUpInitChildren(csmc, skip=[csmc.cliq.index;])                                        # FIXME, post #459 calls?
  # remove all lingering upmessage likelihoods
  oldTags = lsf(csmc.cliqSubFg, tags=[:LIKELIHOODMESSAGE;])
  0 < length(oldTags) ? @warn("stale LIKELIHOODMESSAGE tags present in mustInitUpCliq_StateMachine") : nothing
  oldFcts = oldTags .|> x->getFactor(csmc.cliqSubFg, x)
  # add incoming up messages as priors to subfg
  infocsm(csmc, "8f, mustInitUpCliq_StateMachine -- adding up message factors")
  deleteMsgFactors!(csmc.cliqSubFg, oldFcts)
  # interally adds :LIKELIHOODMESSAGE, :UPWARD_DIFFERENTIAL, :UPWARD_COMMON to each of the factors
  msgfcts = addMsgFactors!(csmc.cliqSubFg, upmsgs, UpwardPass)

  # store the cliqSubFg for later debugging
  opts = getSolverParams(csmc.dfg)
  if opts.dbg
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforeupsolve"))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforeupsolve.pdf"))
  end

  # attempt initialize if necessary
  if !areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq)
    # structure for all up message densities computed during this initialization procedure.
    varorder = getCliqVarInitOrderUp(csmc.tree, csmc.cliq)
    cycleInitByVarOrder!(csmc.cliqSubFg, varorder, logger=csmc.logger)
    # is clique fully upsolved or only partially?
    # print out the partial init status of all vars in clique
    printCliqInitPartialInfo(csmc.cliqSubFg, csmc.cliq, csmc.logger)
  end

  # check again if all cliq vars have been initialized so that full inference can occur on clique
  if areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq)
    infocsm(csmc, "8f, mustInitUpCliq_StateMachine -- all initialized")
    # go to 8g.
    return doCliqUpSolveInitialized_StateMachine
  else
    infocsm(csmc, "8f, mustInitUpCliq_StateMachine -- not able to init all")
    # TODO Simplify this
    status = getCliqueStatus(csmc.cliq)
    status = (status == :initialized || length(getParent(csmc.tree, csmc.cliq)) == 0) ? status : :needdownmsg
    # notify of results (big part of #459 consolidation effort)
    prepPutCliqueStatusMsgUp!(csmc, status)
    # go to 8h
    return rmUpLikeliSaveSubFg_StateMachine
  end
end


"""
    $SIGNATURES

Find upward Chapman-Kolmogorov transit integral solution (approximation).

Notes
- State machine function nr. 8g
- Assumes LIKELIHOODMESSAGE factors are in csmc.cliqSubFg but does not remove them.
- TODO: Make multi-core

DevNotes
- NEEDS DFG v0.8.1, see IIF #760
"""
function doCliqUpSolveInitialized_StateMachine(csmc::CliqStateMachineContainer)

  # check if all cliq vars have been initialized so that full inference can occur on clique
  status = getCliqueStatus(csmc.cliq)
  infocsm(csmc, "8g, doCliqUpSolveInitialized_StateMachine -- clique status = $(status)")
  setCliqDrawColor(csmc.cliq, "red")
  # get Dict{Symbol, TreeBelief} of all updated variables in csmc.cliqSubFg
  retdict = approxCliqMarginalUp!(csmc, logger=csmc.logger)

  # set clique color accordingly, using local memory
  updateFGBT!(csmc.cliqSubFg, csmc.cliq, retdict, dbg=getSolverParams(csmc.cliqSubFg).dbg, logger=csmc.logger) # urt
  setCliqDrawColor(csmc.cliq, isCliqFullDim(csmc.cliqSubFg, csmc.cliq) ? "pink" : "tomato1")

  # notify of results (part of #459 consolidation effort)
  getCliqueData(csmc.cliq).upsolved = true
  status = :upsolved

  # Send upward message, NOTE consolidation WIP #459
  infocsm(csmc, "8g, doCliqUpSolveInitialized_StateMachine -- setting up messages with status = $(status)")
  prepPutCliqueStatusMsgUp!(csmc, status)

  # go to 8h
  return rmUpLikeliSaveSubFg_StateMachine
end


"""
    $SIGNATURES

Close out up solve attempt by removing any LIKELIHOODMESSAGE and save a debug cliqSubFg.

Notes
- State machine function nr. 8h
- Assumes LIKELIHOODMESSAGE factors are in csmc.cliqSubFg and also removes them.
- TODO: Make multi-core

DevNotes
- NEEDS DFG v0.8.1, see IIF #760
"""
function rmUpLikeliSaveSubFg_StateMachine(csmc::CliqStateMachineContainer)
  #
  status = getCliqueStatus(csmc.cliq)

  # remove msg factors that were added to the subfg
  tags__ = getSolverParams(csmc.cliqSubFg).useMsgLikelihoods ? [:UPWARD_COMMON;] : [:LIKELIHOODMESSAGE;]
  msgfcts = lsf(csmc.cliqSubFg, tags=tags__) .|> x->getFactor(csmc.cliqSubFg, x)
  infocsm(csmc, "8g, doCliqUpsSolveInit.! -- status = $(status), removing $(tags__) factors, length=$(length(msgfcts))")
  deleteMsgFactors!(csmc.cliqSubFg, msgfcts)

  # store the cliqSubFg for later debugging
  opts = getSolverParams(csmc.dfg)
  if opts.dbg
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_afterupsolve"))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_afterupsolve.pdf"))
  end

  # go to 9
  return finishCliqSolveCheck_StateMachine
end



"""
$SIGNATURES

WIP to resolve 459 dwnMsg consolidation.  This is partly doing some kind of downsolve but seems out of place.

Notes
- State machine function nr. 10a

DevNotes
- FIXME, resolve/consolidate whatever is going on here
"""
function somebodyLovesMe_StateMachine(csmc::CliqStateMachineContainer)
  infocsm(csmc, "10a, canCliqDownSolve_StateMachine, going to block on parent.")
  prnt = getParent(csmc.tree, csmc.cliq)
  
  # block here until parent is downsolved
  setCliqDrawColor(csmc.cliq, "turquoise")
  # this part is a pull model #674
  while fetchMsgDwnInit(prnt[1]).status != :downsolved
    wait(getSolveCondition(prnt[1]))
  end
  # blockMsgDwnUntilStatus(prnt[1], :downsolved)
  # blockCliqUntilParentDownSolved(, logger=csmc.logger)

  # yes, continue with downsolve
  prntst = getCliqueStatus(prnt[1])
  infocsm(csmc, "10a, somebodyLovesMe_StateMachine, parent status=$prntst.")
  if prntst != :downsolved
    infocsm(csmc, "10a, somebodyLovesMe_StateMachine, going around again.")
    return canCliqDownSolve_StateMachine
  end

  infocsm(csmc, "10a, somebodyLovesMe_StateMachine, going for down solve.")
  # go to 11
  return doCliqDownSolve_StateMachine
end


"""
    $SIGNATURES

Nedd description for this???

Notes
- State machine function nr. 8c

DevNotes
- Must consolidate as part of #459
"""
function waitChangeOnParentCondition_StateMachine(csmc::CliqStateMachineContainer)
  # 
  setCliqDrawColor(csmc.cliq, "coral")

  prnt = getParent(csmc.tree, csmc.cliq)
  if length(prnt) > 0
    infocsm(csmc, "8c, waitChangeOnParentCondition_StateMachine, wait on parent=$(prnt[1].index) for condition notify.")
    @sync begin
      @async begin
        sleep(1)
        notify(getSolveCondition(prnt[1]))
      end
      # wait but don't clear what is in the Condition (guess)
      wait(getSolveCondition(prnt[1]))
    end
  else
    infocsm(csmc, "8c, waitChangeOnParentCondition_StateMachine, cannot wait on parent for condition notify.")
    @warn "no parent!"
  end

  # go to 4b 
  return trafficRedirectConsolidate459_StateMachine
    # # go to 4
    # return canCliqMargSkipUpSolve_StateMachine
end


"""
    $SIGNATURES

WIP #459 dwnMsg consolidation towards blocking cliq that `:needdwninit` to wait on parent `:initialized` dwn message. 

Notes
- State machine function nr.6e

DevNotes
- Seems really unnecessary
- Separated out during #459 dwnMsg consolidation
- Should only happen in long downinit chains below parent that needed dwninit
- TODO figure out whats different between this and 8c
"""
function slowOnPrntAsChildrNeedDwn_StateMachine(csmc::CliqStateMachineContainer)
  # do actual fetch
  prtmsg = fetchMsgDwnInit(getParent(csmc.tree, csmc.cliq)[1]).status
  if prtmsg == :initialized
    # FIXME what does this mean???
    # probably that downward init should commence (not complete final solve)
      allneeddwn = true
    # go to 8e.ii
    # return attemptDownInit_StateMachine
  end

  # FIXME WHY THIS???
  # go to 7
  return determineCliqNeedDownMsg_StateMachine
end


"""
    $SIGNATURES

Decide whether to pursue and upward or downward solve with present state.

Notes
- State machine function nr. 7c
"""
function towardUpOrDwnSolve_StateMachine(csmc::CliqStateMachineContainer)

  sleep(0.1) # FIXME remove after #459 resolved
  # return doCliqInferAttempt_StateMachine
  cliqst = getCliqueStatus(csmc.cliq)
  infocsm(csmc, "7c, status=$(cliqst), before picking direction")
  # d1,d2,cliqst = doCliqInitUpOrDown!(csmc.cliqSubFg, csmc.tree, csmc.cliq, isprntnddw)
  if cliqst == :needdownmsg && !isCliqParentNeedDownMsg(csmc.tree, csmc.cliq, csmc.logger)
      # FIXME, end 459, was collectDwnInitMsgFromParent_StateMachine but 674 requires pull model
      # go to 8c
      # notifyCSMCondition(csmc.cliq)
      # return waitChangeOnParentCondition_StateMachine
    # go to 8a
    return collectDwnInitMsgFromParent_StateMachine
  # HALF DUPLICATED IN STEP 4
  elseif cliqst == :marginalized
    # go to 1
    return isCliqUpSolved_StateMachine
  end

  # go to 8b
  return attemptCliqInitUp_StateMachine
end


"""
    $SIGNATURES

Quick redirection of out-marginalized cliques to downsolve path, or wait on children cliques to get a csm status.

Notes
- State machine function nr.4
"""
function canCliqMargSkipUpSolve_StateMachine(csmc::CliqStateMachineContainer)

  cliqst = getCliqueStatus(csmc.oldcliqdata)
  infocsm(csmc, "4, canCliqMargSkipUpSolve_StateMachine, $cliqst, csmc.incremental=$(csmc.incremental)")

  # if clique is out-marginalized, then no reason to continue with upsolve
  # marginalized state is set in `testCliqCanRecycled_StateMachine`
  if cliqst == :marginalized
    # go to 10 -- Add case for IIF issue #474
    return canCliqDownSolve_StateMachine
  end

  # go to 4e
  return blockUntilChildrenHaveStatus_StateMachine
end


"""
    $SIGNATURES

Do down solve calculations, loosely translates to solving Chapman-Kolmogorov
transit integral in downward direction.

Notes
- State machine function nr. 8e.ii.
  - Follows routines in 8c.
- Pretty major repeat of functionality, FIXME
- TODO: Make multi-core

Initialization requires down message passing of more specialized down init msgs.
This function performs any possible initialization of variables and retriggers
children cliques that have not yet initialized.

Notes:
- Assumed this function is only called after status from child clique up inits completed.
- Assumes cliq has parent.
  - will fetch message from parent
- Will perform down initialization if status == `:needdownmsg`.
- might be necessary to pass furhter down messges to child cliques that also `:needdownmsg`.
- Will not complete cliq solve unless all children are `:upsolved` (upward is priority).
- `dwinmsgs` assumed to come from parent initialization process.
- assume `subfg` as a subgraph that can be modified by this function (add message factors)
  - should remove message prior factors from subgraph before returning.
- May modify `cliq` values.
  - `putMsgUpInit!(cliq, msg)`
  - `setCliqueStatus!(cliq, status)`
  - `setCliqDrawColor(cliq, "sienna")`
  - `notifyCliqDownInitStatus!(cliq, status)`

Algorithm:
- determine which downward messages influence initialization order
- initialize from singletons to most connected non-singletons
- revert back to needdownmsg if cycleInit does nothing
- can only ever return :initialized or :needdownmsg status

DevNotes
- TODO Lots of cleanup required, especially from calling function.
- TODO move directly into a CSM state function
"""
function attemptDownInit_StateMachine(csmc::CliqStateMachineContainer)
  setCliqDrawColor(csmc.cliq, "green")
  opt = getSolverParams(csmc.cliqSubFg)
  dwnkeys_ = lsf(csmc.cliqSubFg, tags=[:DOWNWARD_COMMON;]) .|> x->ls(csmc.cliqSubFg, x)[1]

  ## TODO deal with partial inits only, either delay or continue at end...
  # find intersect between downinitmsgs and local clique variables
  # if only partials available, then
  infocsm(csmc, "8e.ii, attemptDownInit_StateMachine, do cliq init down dwinmsgs=$(dwnkeys_)")
  # get down variable initialization order
  initorder = getCliqInitVarOrderDown(csmc.cliqSubFg, csmc.cliq, dwnkeys_)
  infocsm(csmc, "8e.ii, attemptDownInit_StateMachine,  initorder=$(initorder)")

  # store the cliqSubFg for later debugging
  if opt.dbg
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opt.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforedowninit"))
  end

  # cycle through vars and attempt init
  infocsm(csmc, "8e.ii, attemptDownInit_StateMachine, cycle through vars and attempt init")
  cliqst = :needdownmsg
  if cycleInitByVarOrder!(csmc.cliqSubFg, initorder)
    cliqst = :initialized
  end

  infocsm(csmc, "8e.ii, attemptDownInit_StateMachine, current status: $cliqst")
  # store the cliqSubFg for later debugging
  if opt.dbg
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opt.logpath,"logs/cliq$(csmc.cliq.index)/fg_afterdowninit"))
  end

  # TODO: transfer values changed in the cliques should be transfered to the tree in proc 1 here.
  # # TODO: is status of notify required here?
  setCliqueStatus!(csmc.cliq, cliqst)

  # go to 8l
  return rmMsgLikelihoodsAfterDwn_StateMachine
end

"""
    $SIGNATURES

Remove any `:LIKELIHOODMESSAGE` from `cliqSubFg`.

Notes
- State machine function nr.8l
"""
function rmMsgLikelihoodsAfterDwn_StateMachine(csmc::CliqStateMachineContainer)
  ## TODO only remove :DOWNWARD_COMMON messages here
  #
  ## FIXME move this to separate state in CSM.
  # remove all message factors
  # remove msg factors previously added
  fctstorm = ls(csmc.cliqSubFg, tags=[:LIKELIHOODMESSAGE;])
  infocsm(csmc, "8e.ii., attemptDownInit_StateMachine, removing factors $fctstorm")
  rmfcts = fctstorm .|> x->getFactor(csmc.cliqSubFg, x)
  deleteMsgFactors!(csmc.cliqSubFg, rmfcts )

  # go to 8d
  return decideUpMsgOrInit_StateMachine
end

"""
    $SIGNATURES

Blocking until all children of the clique's parent (i.e. siblings) have a valid status.

Notes
- State machine function nr. 5
"""
function blockUntilSiblingsStatus_StateMachine(csmc::CliqStateMachineContainer)
  infocsm(csmc, "5, blocking on parent until all sibling cliques have valid status")
  setCliqDrawColor(csmc.cliq, "blueviolet")

  cliqst = getCliqueStatus(csmc.cliq)
  infocsm(csmc, "5, block on siblings")
  prnt = getParent(csmc.tree, csmc.cliq)
  if length(prnt) > 0
    infocsm(csmc, "5, has parent clique=$(prnt[1].index)")
    ret = fetchChildrenStatusUp(csmc.tree, prnt[1], csmc.logger)
    infocsm(csmc,"prnt $(prnt[1].index), fetched all, keys=$(keys(ret)).")
  end

  infocsm(csmc, "5, finishing")
  # go to 6c
  return doesParentNeedDwn_StateMachine
end


"""
    $SIGNATURES

Delay loop if waiting on upsolves to complete.

Notes
- State machine 7b
- Also has "recursion", to come back to this function to make sure that child clique updates are in fact upsolved.
- Differs from 4e in that here children must be "upsolved" or equivalent to continue.
"""
function slowIfChildrenNotUpSolved_StateMachine(csmc::CliqStateMachineContainer)

  # childs = getChildren(csmc.tree, csmc.cliq)
  # len = length(childs)
  @inbounds for chld in getChildren(csmc.tree, csmc.cliq)
    chst = getCliqueStatus(chld)
    if !(chst in [:upsolved;:uprecycled;:marginalized])
      infocsm(csmc, "7b, slowIfChildrenNotUpSolved_StateMachine, wait condition on upsolve, cliq=$(chld.index), ch_lbl=$(getCliqFrontalVarIds(chld)[1]).")
      
      # wait for child clique status/msg to be updated
      wait(getSolveCondition(chld))
      # tsk = @async 
      # timedwait(()->tsk.state==:done, 10)

      # check again and reroute if :needdownmsg
      # chst = getCliqueStatus(chld)
      # if chst == :needdownmsg
      #   return prntPrepDwnInitMsg_StateMachine
      # end
    end
  end

  # go to 4b
  return trafficRedirectConsolidate459_StateMachine
end


# """
#     $SIGNATURES

# New function in 459 dwnMsg init for when child clique needdowninit msg.
  
# Notes
# - State machine function nr. 7d

# DevNotes
# - "Parent" must do work from (pre-#459) `collectDwnInitMsgFromParent_StateMachine`.
# """
# function prntPrepDwnInitMsg_StateMachine(csmc::CliqStateMachineContainer)
  
#   @warn("prntPrepDwnInitMsg_StateMachine -- WIP")
#   infocsm(csmc, "prntPrepDwnInitMsg_StateMachine -- WIP")

#   # TEMP go to 7b
#   return slowIfChildrenNotUpSolved_StateMachine
# end



























"""
    $SIGNATURES

Block until all children have a csm status, using `fetch---Condition` model. 

Notes
- State machine function nr.4e
"""
function blockUntilChildrenHaveStatus_StateMachine(csmc::CliqStateMachineContainer)
  #must happen before if :null
  stdict = fetchChildrenStatusUp(csmc.tree, csmc.cliq, csmc.logger)
  infocsm(csmc,"fetched all, keys=$(keys(stdict)).")
  # make sure forceproceed is false, not strictly needed if csmc starts false
  csmc.forceproceed = false # only used in 7

  # go to 4b
  return trafficRedirectConsolidate459_StateMachine
end



## ============================================================================================
# start of things downinit
## ============================================================================================


"""
$SIGNATURES

Test waiting order between siblings for cascading downward tree initialization.
  
Notes
- State machine function 8j.

DevNotes
- FIXME, something wrong with CSM sequencing, https://github.com/JuliaRobotics/IncrementalInference.jl/issues/602#issuecomment-682114232
- This might be replaced with 4-stroke tree-init, if that algorithm turns out to work in all cases.
"""
function dwnInitSiblingWaitOrder_StateMachine(csmc::CliqStateMachineContainer)
  
  prnt = getParent(csmc.tree, csmc.cliq)[1]
  opt = getSolverParams(csmc.cliqSubFg) # csmc.dfg
  
  # now get the newly computed message from the appropriate container
  # make sure this is a pull model #674 (pull msg from source/prnt)
  # FIXME must be consolidated as part of #459
  dwinmsgs = getfetchCliqueInitMsgDown(getCliqueData(prnt), from=:dwnInitSiblingWaitOrder_StateMachine)

  # add downward belief prop msgs
  msgfcts = addMsgFactors!(csmc.cliqSubFg, dwinmsgs, DownwardPass)
  # determine if more info is needed for partial
  sdims = getCliqVariableMoreInitDims(csmc.cliqSubFg, csmc.cliq)
  updateCliqSolvableDims!(csmc.cliq, sdims, csmc.logger)
    
  opt.dbg ? saveDFG(joinLogPath(csmc.cliqSubFg, "logs", "cliq$(csmc.cliq.index)", "fg_DWNCMN"), csmc.cliqSubFg) : nothing

  dwnkeys_ = collect(keys(dwinmsgs.belief))

  # NOTE, only use separators, not all parent variables
  # dwnkeys_ = lsf(csmc.cliqSubFg, tags=[:DOWNWARD_COMMON;]) .|> x->ls(csmc.cliqSubFg, x)[1]
  # @assert length(intersect(dwnkeys, dwnkeys_)) == length(dwnkeys) "split dwnkeys_ is not the same, $dwnkeys, and $dwnkeys_"

  # priorize solve order for mustinitdown with lowest dependency first
  # follow example from issue #344
  mustwait = false
  if length(intersect(dwnkeys_, getCliqSeparatorVarIds(csmc.cliq))) == 0
    infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, no can do, must wait for siblings to update parent first.")
    mustwait = true
  elseif getSiblingsDelayOrder(csmc.tree, csmc.cliq, dwnkeys_, logger=csmc.logger)
    infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, prioritize")
    mustwait = true
  elseif getCliqSiblingsPartialNeeds(csmc.tree, csmc.cliq, dwinmsgs, logger=csmc.logger)
    infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, partialneedsmore")
    mustwait = true
  end

  solord = getCliqSiblingsPriorityInitOrder( csmc.tree, prnt, csmc.logger )
  noOneElse = areSiblingsRemaingNeedDownOnly(csmc.tree, csmc.cliq)
  infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, $(prnt.index), $mustwait, $noOneElse, solord = $solord")

  if mustwait && csmc.cliq.index != solord[1] # && !noOneElse
    infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, must wait on change.")
    # remove all message factors
    fctstorm = ls(csmc.cliqSubFg, tags=[:DOWNWARD_COMMON;])
    infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, removing factors $fctstorm")
    rmfcts = fctstorm .|> x->getFactor(csmc.cliqSubFg, x)
    # remove msg factors previously added
    deleteMsgFactors!(csmc.cliqSubFg, rmfcts )

    # go to 8c
    return waitChangeOnParentCondition_StateMachine
  end

  # go to 8e.ii.
  return attemptDownInit_StateMachine
end


"""
    $SIGNATURES

Do down initialization calculations, loosely translates to solving Chapman-Kolmogorov
transit integral in downward direction.

Notes
- State machine function nr. 8a
- Includes initialization routines.
- TODO: Make multi-core

DevNotes
- FIXME major refactor of this function required.
- FIXME this function actually occur during the parent CSM, therefore not all pull model #674
"""
function collectDwnInitMsgFromParent_StateMachine(csmc::CliqStateMachineContainer)
  #
  # TODO consider exit early for root clique rather than avoiding this function
  infocsm(csmc, "8a, needs down message -- attempt down init")
  setCliqDrawColor(csmc.cliq, "gold")

  # initialize clique in downward direction
  # not if parent also needs downward init message
  # prnt = getParent(csmc.tree, csmc.cliq)[1]
  opt = getSolverParams(csmc.dfg)
  @assert !haskey(opt.devParams,:dontUseParentFactorsInitDown) "dbgnew is old school, 459 dwninit consolidation has removed the option for :dontUseParentFactorsInitDown"

  # take atomic lock OF PARENT ??? when waiting for downward information
  # lockUpStatus!(prnt, prnt.index, true, csmc.logger, true, "cliq$(csmc.cliq.index)") # TODO XY ????
  # infocsm(csmc, "8a, after up lock")

  # get down message from the parent
  # check if any msgs should be multiplied together for the same variable
  # get the current messages ~~stored in~~ [going to] the parent (pull model #674)
  # FIXME, post #459 calls?
  # this guy is getting any sibling up messages by calling on the parent
  prntmsgs::Dict{Int, LikelihoodMessage} = getMsgsUpInitChildren(csmc.tree, csmc.cliq, TreeBelief, skip=[csmc.cliq.index;])         
  
  # reference to default dict location
  dwinmsgs = getfetchCliqueInitMsgDown(csmc.cliq.data, from=:getMsgDwnThisInit) |> deepcopy  #JT 459 products = getMsgDwnThisInit(prnt)
  infocsm(csmc, "getMsgInitDwnParent -- msg ids::Int=$(collect(keys(prntmsgs)))")
  
  # stack all parent incoming upward messages into dict of vector msgs
  prntBelDictVec::Dict{Symbol, Vector{TreeBelief}} = convertLikelihoodToVector(prntmsgs, logger=csmc.logger)
  ## TODO use parent factors too
  # intersect with the asking clique's separator variables
  # this function populates `dwinmsgs` with the appropriate products described in `prntBelDictVec`
  # FIXME, should not be using full .dfg ???
  condenseDownMsgsProductPrntFactors!(csmc.dfg, dwinmsgs, prntBelDictVec, prnt, csmc.cliq, csmc.logger)

  # remove msgs that have no data
  rmlist = Symbol[]
  for (prsym,belmsg) in dwinmsgs.belief
    if belmsg.inferdim < 1e-10
      # no information so remove
      push!(rmlist, prsym)
    end
  end
  infocsm(csmc, "prepCliqInitMsgsDown! -- rmlist, no inferdim, keys=$(rmlist)")
  for pr in rmlist
    delete!(dwinmsgs.belief, pr)
  end

  infocsm(csmc, "prepCliqInitMsgsDown! -- product keys=$(collect(keys(dwinmsgs.belief)))")

  # now put the newly computed message in the appropriate container
  # FIXME THIS IS A PUSH MODEL, see #674 -- must make pull model first
  # FIXME must be consolidated as part of #459
  putCliqueInitMsgDown!(getCliqueData(csmc.cliq), dwinmsgs)
  # unlock
  # unlockUpStatus!(prnt) # TODO XY ????, maybe can remove after pull model #674?
  infocsm(csmc, "8a, attemptCliqInitD., unlocked")

  # go to 7b (maybe, and part of dwnMsg #459 WIP 9)
  return slowIfChildrenNotUpSolved_StateMachine
    # # go to 8j.
    # return dwnInitSiblingWaitOrder_StateMachine
end


















"""
    $SIGNATURES

Redirect CSM traffic in various directions 

Notes
- State machine function nr.4b
- Was refactored during #459 dwnMsg effort.

DevNotes
- Consolidate with 7?
"""
function trafficRedirectConsolidate459_StateMachine(csmc::CliqStateMachineContainer)

  cliqst = getCliqueStatus(csmc.cliq)
  infocsm(csmc, "4b, trafficRedirectConsolidate459_StateMachine, cliqst=$cliqst")

  # if no parent or parent will not update

  # for recycle computed clique values case
  if csmc.incremental && cliqst == :downsolved
    csmc.incremental = false
    # might be able to recycle the previous clique solve, go to 0b
    return checkChildrenAllUpRecycled_StateMachine
  end

  prnt = getParent(csmc.tree, csmc.cliq)
  if 0 == length(prnt) || cliqst == :needdownmsg
    # go to 7
    return determineCliqNeedDownMsg_StateMachine
  end

  # Some traffic direction
  if cliqst == :null
    # go to 4d
    return checkIfCliqNullBlock_StateMachine
  # elseif cliqst == :needdownmsg  # && cliqst != :null
  #   # got to 4c (seems like only needdownmsg case gets here)
  #   return doAnyChildrenNeedDwn_StateMachine
  else
    # go to 6c
    return doesParentNeedDwn_StateMachine
  end
end
# :null => checkIf
# :needd => untilDown
# :initi => blockCliq
# :upsol => blockCliq




## ============================================================================================
# START of dwnmsg consolidation bonanza
## ============================================================================================
#
# blockCliqSiblingsParentChildrenNeedDown_  # Blocking case when all siblings and parent :needdownmsg.
# doAllSiblingsNeedDwn_StateMachine   # Trying to figure out when to block on siblings for cascade down init. 
# checkIfCliqNullBlock_StateMachine         # If all children (then also escalate to) :needdownmsgs and block until sibling status.
# determineCliqNeedDownMsg_StateMachine     # Try decide whether this `csmc.cliq` needs a downward initialization message.
# doAnyChildrenNeedDwn_StateMachine         # Determine if any one of the children :needdownmsg.
# downInitRequirement_StateMachine          # Place fake up msg and notify down init status if any children :needdownmsg


"""
$SIGNATURES

Trying to figure out when to block on siblings for cascade down init.  

Notes
- State machine function nr.6d
- Part of #459 dwnMsg consolidation work.
- used for regulating long need down message chains.
- exit strategy is parent becomes status `:initialized`.
- Assume there must be a parent.

DevNotes
- Consolidation with work with similar likely required.
"""
function doAllSiblingsNeedDwn_StateMachine(csmc::CliqStateMachineContainer)
    
  prnt = getParent(csmc.tree, csmc.cliq)
  prnt_ = prnt[1]
  
  allneeddwn = true
  for ch in getChildren(csmc.tree, prnt_)
    chst = getCliqueStatus(ch)
    if chst != :needdownmsg
      allneeddwn = false
      break;
    end
  end

  # FIXME, understand why is there another status event from parent msg here... How to consolidate this with CSM 8a
  if allneeddwn
    # go to 6e
    return slowOnPrntAsChildrNeedDwn_StateMachine
  end

  # go to 7
  return determineCliqNeedDownMsg_StateMachine
end





"""
$SIGNATURES

If all children (then also escalate to) :needdownmsgs and block until sibling status.

Notes
- State machine function nr.4d

DevNotes
- TODO Any overlap with nr.4c??
"""
function checkIfCliqNullBlock_StateMachine(csmc::CliqStateMachineContainer)
  # fetch (should not block)
  stdict = fetchChildrenStatusUp(csmc.tree, csmc.cliq, csmc.logger)
  infocsm(csmc,"fetched all, keys=$(keys(stdict)).")
  chstatus = collect(values(stdict))
  len = length(chstatus)

  # if all children needdownmsg
  if len > 0 && sum(chstatus .== :needdownmsg) == len
    # TODO maybe can happen where some children need more information?
    infocsm(csmc, "4d, checkIfCliqNullBlock_StateMachine, escalating to :needdownmsg since all children :needdownmsg")

    # NOTE, trying consolidation with prepPutUp for #459 effort
    prepPutCliqueStatusMsgUp!(csmc, :needdownmsg)
    setCliqDrawColor(csmc.cliq, "yellowgreen")

    # debuggin #459 transition
    infocsm(csmc, "4d, checkIfCliqNullBlock_StateMachine -- finishing before going to  blockUntilSiblingsStatus_StateMachine")

    # go to 5
    return blockUntilSiblingsStatus_StateMachine
  end

  # go to 6c
  return doesParentNeedDwn_StateMachine
end


"""
    $SIGNATURES

Try decide whether this `csmc.cliq` needs a downward initialization message.

Notes
- State machine function nr. 7

DevNotes
- Consolidate with 4b?
"""
function determineCliqNeedDownMsg_StateMachine(csmc::CliqStateMachineContainer)

  # # # old 4c
  # # function doAnyChildrenNeedDwn_StateMachine(csmc::CliqStateMachineContainer)
  #   for ch in getChildren(csmc.tree, csmc.cliq)
  #     if getCliqueStatus(ch) == :needdownmsg
  #       infocsm(csmc, "4c, doAnyChildrenNeedDwn_StateMachine, must deal with child :needdownmsg")
  #       csmc.forceproceed = true
  #       break
  #     end
  #   end
  
  #   if !csmc.forceproceed
  #     infocsm(csmc, "4c, doAnyChildrenNeedDwn_StateMachine, no children :needdownmsg")
  #     # go to 5
  #     return blockUntilSiblingsStatus_StateMachine
  #   end
  #   infocsm(csmc, "4c, doAnyChildrenNeedDwn_StateMachine, yes some children do :needdownmsg")
  
  #   # go to 6c
  #   return doesParentNeedDwn_StateMachine
  # # end

  # FIXME add path to 8a
  # return collectDwnInitMsgFromParent_StateMachine



  # fetch children status
  stdict = fetchChildrenStatusUp(csmc.tree, csmc.cliq, csmc.logger)
  infocsm(csmc,"fetched all, keys=$(keys(stdict)).")
  
  # hard assumption here on upsolve from leaves to root
  hasinit = false
  resolveinit = true
  # fetch status from children (should already be available -- i.e. should not block)
  for (clid, clst) in stdict
    infocsm(csmc, "7, check stdict children: clid=$(clid), clst=$(clst)")
    ## :needdownmsg # 'send' downward init msg direction
    (clst in [:upsolved;:marginalized;:downsolved;:uprecycled]) ? nothing : (resolveinit = false)
    hasinit = clst == :needdownmsg ? true : hasinit
  end
  infocsm(csmc, "7, resolveinit=$(resolveinit)")
  infocsm(csmc, "7, start, forceproceed=$(csmc.forceproceed)")
  
  # merged in from 4c here into 7, part of dwnMsg #459
  if hasinit
    # go to 5
    return blockUntilSiblingsStatus_StateMachine
  end

  if csmc.forceproceed || resolveinit # proceed
    # TODO, remove csmc.forceproceed, only used in 7
    csmc.forceproceed = false
    
    # go to 8j (dwnMsg #459 WIP 9)
    return dwnInitSiblingWaitOrder_StateMachine
      # # go to 7c
      # return towardUpOrDwnSolve_StateMachine
  end

  # special case short cut
  cliqst = getCliqueStatus(csmc.cliq)
  if cliqst == :needdownmsg
    infocsm(csmc, "7b, shortcut since this cliq :needdownmsg.")
    # go to 4b
    return trafficRedirectConsolidate459_StateMachine # canCliqMargSkipUpSolve_StateMachine
  end

  # go to 7b
  return slowIfChildrenNotUpSolved_StateMachine
end


# """
#     $SIGNATURES

# Direct traffic if any one of the children :needdownmsg.

# Notes
# - State machine function nr.4c

# DevNotes
# - Consolidate into 7 
#   - Consolidate with 8d???
# - FIXME, likely location to divert to wait on parent completion of combo dwninitmsg.
#   - in place of `forceproceed`
# - TODO remove csmc.forceproceed entirely from CSM
#   - Try force parent to initialize ?? FIXME DEPRECATE
# """
# function doAnyChildrenNeedDwn_StateMachine(csmc::CliqStateMachineContainer)
#   for ch in getChildren(csmc.tree, csmc.cliq)
#     if getCliqueStatus(ch) == :needdownmsg
#       infocsm(csmc, "4c, doAnyChildrenNeedDwn_StateMachine, must deal with child :needdownmsg")
#       csmc.forceproceed = true
#       break
#     end
#   end

#   if !csmc.forceproceed
#     infocsm(csmc, "4c, doAnyChildrenNeedDwn_StateMachine, no children :needdownmsg")
#     # go to 5
#     return blockUntilSiblingsStatus_StateMachine
#   end
#   infocsm(csmc, "4c, doAnyChildrenNeedDwn_StateMachine, yes some children do :needdownmsg")

#   # go to 6c
#   return doesParentNeedDwn_StateMachine
# end


"""
$SIGNATURES

Place up msg and notify down init status if any children :needdownmsg.

Notes
- StateMachine function nr. 8d
- TODO figure out if this function is duplication of other needdwnmsg functionality?
  - Consolidate with 7???
  - Consolidate with 4c???
"""
function decideUpMsgOrInit_StateMachine(csmc::CliqStateMachineContainer)
  #
  infocsm(csmc, "8d, downInitRequirement_StateMachine., start")
  
  # children = getChildren(csmc.tree, csmc.cliq)
  # if doAnyChildrenNeedDwnMsg(children)
  someChildrenNeedDwn = false
  for ch in getChildren(csmc.tree, csmc.cliq)
    if getCliqueStatus(ch) == :needdownmsg
      someChildrenNeedDwn = true
      break
    end
  end

  if someChildrenNeedDwn
    # go to 8k
    return sendCurrentUpMsg_StateMachine
  end
  
  # go to 8b
  return attemptCliqInitUp_StateMachine
end


"""
    $SIGNATURES

Blocking case when all siblings and parent :needdownmsg.

Notes
- State machine function nr. 6c

DevNotes
- FIXME if statements can be slightly simplified before further consolidation.
- FIXME understand if this should be consolidated with 4b. `trafficRedirectConsolidate459_StateMachine`?
"""
function doesParentNeedDwn_StateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "6c, check/block sibl&prnt :needdownmsg")
  prnt = getParent(csmc.tree, csmc.cliq)
  if 0 == length(prnt) || getCliqueStatus(prnt[1]) != :needdownmsg
    # go to 7
    return determineCliqNeedDownMsg_StateMachine
  # elseif
  end

  # if getCliqueStatus(prnt[1]) != :needdownmsg
  #   # go to 7
  #   return determineCliqNeedDownMsg_StateMachine
  # end
  
  # go to 6d
  return doAllSiblingsNeedDwn_StateMachine
end


"""
    $SIGNATURES

Determine if up initialization calculations should be attempted.

Notes
- State machine function nr. 8b
"""
function attemptCliqInitUp_StateMachine(csmc::CliqStateMachineContainer)
  
  # should calculations be avoided.
  notChildNeedDwn = !doAnyChildrenNeedDwnMsg(csmc.tree, csmc.cliq) 
  infocsm(csmc, "8b, attemptCliqInitUp, !doAnyChildrenNeedDwnMsg()=$(notChildNeedDwn)" )
  if getCliqueStatus(csmc.cliq) in [:initialized; :null; :needdownmsg] && notChildNeedDwn
    # go to 8f.
    return mustInitUpCliq_StateMachine
  end

  # go to 9
  return finishCliqSolveCheck_StateMachine
end


## ============================================================================================
# End of dwnmsg consolidation bonanza
## ============================================================================================


"""
    $SIGNATURES

Placeholder function part of #459 dwnMsg consolidation to send current up message, part of :needdownmsg downinit cascading.

Notes
- State machine function nr. 8k
"""
function sendCurrentUpMsg_StateMachine(csmc::CliqStateMachineContainer)
  # set messages if children :needdownmsg
  infocsm(csmc, "8k, sendCurrentUpMsg_StateMachine -- must set messages for future down init")
  # construct init's up msg to place in parent from initialized separator variables

  # consolidated up messaging (#459)
  infocsm(csmc, "8k, sendCurrentUpMsg_StateMachine -- putting fake upinitmsg in this cliq")
  upmsg = prepCliqInitMsgsUp(csmc.cliqSubFg, csmc.cliq, getCliqueStatus(csmc.cliq))
  prepPutCliqueStatusMsgUp!(csmc, upmsg=upmsg)

  # also send a down message -- seem weird while doing #459 but okay
  cliqst = getCliqueStatus(csmc.cliq)
  notifyCliqDownInitStatus!(csmc.cliq, cliqst, logger=csmc.logger)

  # Legend: initialized but not solved yet (likely child cliques that depend on downward autoinit msgs),
  setCliqDrawColor(csmc.cliq, "sienna")

  infocsm(csmc, "8k, sendCurrentUpMsg_StateMachine -- near-end down init attempt, $cliqst.")

  # go to 8b
  return attemptCliqInitUp_StateMachine
end







































"""
    $SIGNATURES

Build a sub factor graph for clique variables from the larger factor graph.

Notes
- State machine function nr.2
"""
function buildCliqSubgraph_StateMachine(csmc::CliqStateMachineContainer)
  # build a local subgraph for inference operations
  infocsm(csmc, "2, build subgraph syms=$(getCliqAllVarIds(csmc.cliq))")
  buildCliqSubgraph!(csmc.cliqSubFg, csmc.dfg, csmc.cliq)

  # if dfg, store the cliqSubFg for later debugging
  dbgSaveDFG(csmc.cliqSubFg, "cliq$(csmc.cliq.index)/fg_build")

  # go to 4
  return canCliqMargSkipUpSolve_StateMachine
end

"""
    $SIGNATURES

Build a sub factor graph for clique variables from the larger factor graph.

Notes
- State machine function nr.2r
"""
function buildCliqSubgraphForDown_StateMachine(csmc::CliqStateMachineContainer)
  # build a local subgraph for inference operations
  syms = getCliqAllVarIds(csmc.cliq)
  infocsm(csmc, "2r, build subgraph syms=$(syms)")
  csmc.cliqSubFg = buildSubgraph(csmc.dfg, syms, 1; verbose=false)

  opts = getSolverParams(csmc.dfg)
  # store the cliqSubFg for later debugging
  if opts.dbg
    mkpath(joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)"))
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_build_down"))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_build_down.pdf"))
  end

  # go to 10
  return canCliqDownSolve_StateMachine
end

"""
    $SIGNATURES

Either construct and notify of a new upward initialization message and progress to downsolve checks,
or circle back and start building the local clique subgraph.

Notes
- State machine function nr.1
- Root clique message should be empty since it has an empty separator.
"""
function isCliqUpSolved_StateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "1, isCliqUpSolved_StateMachine")
  cliqst = getCliqueStatus(csmc.cliq)

  # if upward complete for any reason, prepare and send new upward message
  if cliqst in [:upsolved; :downsolved; :marginalized; :uprecycled]
    # construct init's up msg from initialized separator variables
    # NOTE cliqSubFg has not been copied yet
    prepPutCliqueStatusMsgUp!(csmc, cliqst, dfg=csmc.dfg)
    #go to 10
    return canCliqDownSolve_StateMachine
  end
  # go to 2
  return buildCliqSubgraph_StateMachine
end


"""
    $SIGNATURES

Final determination on whether can promote clique to `:uprecycled`.

Notes
- State machine function nr.0b
- Assume children clique status is available
- Will return to regular init-solve if new information in children -- ie not uprecycle or marginalized
"""
function checkChildrenAllUpRecycled_StateMachine(csmc::CliqStateMachineContainer)
  count = Int[]
  chldr = getChildren(csmc.tree, csmc.cliq)
  for ch in chldr
    chst = getCliqueStatus(ch)
    if chst in [:uprecycled; :marginalized]
      push!(count, 1)
    end
  end
  infocsm(csmc, "0b, checkChildrenAllUpRecycled_StateMachine -- length(chldr)=$(length(chldr)), sum(count)=$(sum(count))")

  # all children can be used for uprecycled -- i.e. no children have new information
  if sum(count) == length(chldr)
    # set up msg and exit go to 1
    sdims = Dict{Symbol,Float64}()
    for varid in getCliqAllVarIds(csmc.cliq)
      sdims[varid] = 0.0
    end

    # NOTE busy consolidating #459
    updateCliqSolvableDims!(csmc.cliq, sdims, csmc.logger)
        # setCliqueStatus!(csmc.cliq, :uprecycled)
        # replacing similar functionality from CSM 1.
    if getSolverParams(csmc.cliqSubFg).dbg
      tmnow = now()
      tmpst = getCliqueStatus(csmc.cliq)
      @async begin
        mkpath(joinLogPath(csmc.cliqSubFg,"logs","cliq$(csmc.cliq.index)"))
        open(joinLogPath(csmc.cliqSubFg,"logs","cliq$(csmc.cliq.index)","incremental.log"), "w") do f
          println(f, tmnow, ", marginalized from previous status ", tmpst)
        end
      end
    end
    prepPutCliqueStatusMsgUp!(csmc, :uprecycled, dfg=csmc.dfg)
    setCliqDrawColor(csmc.cliq, "orange")
    #go to 10
    return canCliqDownSolve_StateMachine
        # # go to 1
        # return isCliqUpSolved_StateMachine
  end

  # return to regular solve, go to 2
  return buildCliqSubgraph_StateMachine
end

"""
    $SIGNATURES

Determine if clique is upsolved by incremental update and exit the state machine.

Notes
- State machine function nr.0c
- can recycle if two checks:
  - previous clique was identically downsolved
  - all children are also :uprecycled
"""
function testCliqCanIncremtUpdate_StateMachine!(csmc::CliqStateMachineContainer)
	# check if should be trying and can recycle clique computations
    if csmc.incremental && getCliqueStatus(csmc.oldcliqdata) == :downsolved
      csmc.cliq.data.isCliqReused = true
      # check if a subgraph will be needed later
      if csmc.dodownsolve
        # yes need subgraph and need more checks, so go to 2
        return buildCliqSubgraph_StateMachine
      else
        # one or two checks say yes, so go to 4
        return canCliqMargSkipUpSolve_StateMachine
      end
    end

    # nope, regular clique init-solve, go to 1
    return isCliqUpSolved_StateMachine
end

"""
    $SIGNATURES

Notify possible parent if clique is upsolved and exit the state machine.

Notes
- State machine function nr.0
- can recycle if two checks:
  - previous clique was identically downsolved
  - all children are also :uprecycled
"""
function testCliqCanRecycled_StateMachine(csmc::CliqStateMachineContainer)
  # @show getCliqFrontalVarIds(csmc.oldcliqdata), getCliqueStatus(csmc.oldcliqdata)
  infocsm(csmc, "0., $(csmc.incremental) ? :uprecycled => getCliqueStatus(csmc.oldcliqdata)=$(getCliqueStatus(csmc.oldcliqdata))")

  if areCliqVariablesAllMarginalized(csmc.dfg, csmc.cliq)

    # no work required other than assembling upward message
    if getSolverParams(csmc.cliqSubFg).dbg
      tmnow = now()
      tmpst = getCliqueStatus(csmc.cliq)
      @async begin
        mkpath(joinLogPath(csmc.cliqSubFg,"logs","cliq$(csmc.cliq.index)"))
        open(joinLogPath(csmc.cliqSubFg,"logs","cliq$(csmc.cliq.index)","marginalization.log"), "w") do f
          println(f, tmnow, ", marginalized from previous status ", tmpst)
        end
      end
    end
    prepPutCliqueStatusMsgUp!(csmc, :marginalized, dfg=csmc.dfg)

    # set marginalized color
    setCliqDrawColor(csmc.cliq, "blue")

    # set flag, looks to be previously unused???
    getCliqueData(csmc.cliq).allmarginalized = true

    # FIXME divert to rapid CSM exit
	# GUESSING THIS THE RIGHT WAY go to 4
	# return canCliqMargSkipUpSolve_StateMachine
  end

  # go to 0c.
  return testCliqCanIncremtUpdate_StateMachine!
end



## ==============================================================================================
# Launch the FSM
## ==============================================================================================

"""
    $SIGNATURES

Perform upward inference using a state machine solution approach.

Notes:
- will call on values from children or parent cliques
- can be called multiple times
- Assumes all cliques in tree are being solved simultaneously and in similar manner.
- State machine rev.1 -- copied from first TreeBasedInitialization.jl.
- Doesn't do partial initialized state properly yet.
"""
function cliqInitSolveUpByStateMachine!(dfg::G,
                                        tree::AbstractBayesTree,
                                        cliq::TreeClique,
                                        timeout::Union{Nothing, <:Real}=nothing;
                                        N::Int=100,
                                        verbose::Bool=false,
                                        verbosefid=stdout,
                                        oldcliqdata::BayesTreeNodeData=BayesTreeNodeData(),
                                        drawtree::Bool=false,
                                        show::Bool=false,
                                        incremental::Bool=true,
                                        limititers::Int=-1,
                                        upsolve::Bool=true,
                                        downsolve::Bool=true,
                                        recordhistory::Bool=false,
                                        delay::Bool=false,
                                        injectDelayBefore::Union{Nothing,Pair{<:Function, <:Real}}=nothing,
                                        logger::SimpleLogger=SimpleLogger(Base.stdout)) where {G <: AbstractDFG, AL <: AbstractLogger}
  #
  children = getChildren(tree, cliq)#Graphs.out_neighbors(cliq, tree.bt)

  prnt = getParent(tree, cliq)

  destType = (G <: InMemoryDFGTypes) ? G : InMemDFGType

  csmc = CliqStateMachineContainer(dfg, initfg(destType, solverParams=getSolverParams(dfg)), tree, cliq, prnt, children, false, incremental, drawtree, downsolve, delay, getSolverParams(dfg), Dict{Symbol,String}(), oldcliqdata, logger)

  nxt = upsolve ? testCliqCanRecycled_StateMachine : (downsolve ? testCliqCanRecycled_StateMachine : error("must attempt either up or down solve"))

  csmiter_cb = getSolverParams(dfg).drawCSMIters ? ((st::StateMachine)->(cliq.attributes["xlabel"] = st.iter)) : ((st)->())

  statemachine = StateMachine{CliqStateMachineContainer}(next=nxt, name="cliq$(cliq.index)")
  while statemachine(csmc, timeout, verbose=verbose, verbosefid=verbosefid, iterlimit=limititers, recordhistory=recordhistory, housekeeping_cb=csmiter_cb, injectDelayBefore=injectDelayBefore); end
  statemachine.history
end



#
