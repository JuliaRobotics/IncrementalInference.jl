# clique state machine for tree based initialization and inference

# newer exports
# export towardUpOrDwnSolve_StateMachine, maybeNeedDwnMsg_StateMachine, doAnyChildrenNeedDwn_StateMachine
# export prepInitUp_StateMachine, doCliqUpSolveInitialized_StateMachine
# export rmUpLikeliSaveSubFg_StateMachine
# export blockCliqSiblingsParentChildrenNeedDown_StateMachine

export  doCliqDownSolve_StateMachine,
        cleanupAfterDownSolve_StateMachine,
        specialCaseRootDownSolve_StateMachine,
        canCliqDownSolve_StateMachine,
        checkUpsolveFinished_StateMachine,
        prepInitUp_StateMachine,
        doCliqUpSolveInitialized_StateMachine,
        rmUpLikeliSaveSubFg_StateMachine,
        wipRedirect459Dwn_StateMachine,
        waitChangeOnParentCondition_StateMachine,
        slowOnPrntAsChildrNeedDwn_StateMachine,
        towardUpOrDwnSolve_StateMachine,
        canCliqMargSkipUpSolve_StateMachine,
        attemptDownInit_StateMachine,
        rmMsgLikelihoodsAfterDwn_StateMachine,
        blockSiblingStatus_StateMachine,
        slowIfChildrenNotUpSolved_StateMachine,
        blockUntilChildrenHaveStatus_StateMachine,
        dwnInitSiblingWaitOrder_StateMachine,
        trafficRedirectConsolidate459_StateMachine,
        doAllSiblingsNeedDwn_StateMachine,
        maybeNeedDwnMsg_StateMachine,
        determineCliqNeedDownMsg_StateMachine,
        tryInitCliq_StateMachine,
        slowWhileInit_StateMachine,
        doAnyChildrenNeedDwn_StateMachine,
        decideUpMsgOrInit_StateMachine,
        attemptCliqInitUp_StateMachine,
        sendCurrentUpMsg_StateMachine,
        buildCliqSubgraph_StateMachine,
        buildCliqSubgraphForDown_StateMachine,
        isCliqUpSolved_StateMachine,
        checkChildrenAllUpRecycled_StateMachine,
        canCliqIncrRecycle_StateMachine,
        canCliqMargRecycle_StateMachine



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

    prepPutCliqueStatusMsgDwn!(csmc, :downsolved)
    #JT 459 putMsgDwnThis!(cliq, newDwnMsgs), DF still looks like a pull model here #674
    # putMsgDwnThis!(csmc.cliq.data, newDwnMsgs, from=:putMsgDwnThis!) # putCliqueMsgDown!

    # update clique subgraph with new status
    setCliqDrawColor(csmc.cliq, "lightblue")

  csmc.dodownsolve = false
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- finished with downGibbsCliqueDensity, now update csmc")

  # go to 11b.
  return cleanupAfterDownSolve_StateMachine
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
  return wipRedirect459Dwn_StateMachine

end


"""
    $SIGNATURES

Is upsolve complete or should the CSM solving process be repeated.

Notes
- State machine function nr.9

DevNotes
- FIXME FIXME FIXME ensure init worked
"""
function checkUpsolveFinished_StateMachine(csmc::CliqStateMachineContainer)
  cliqst = getCliqueStatus(csmc.cliq)
  infocsm(csmc, "9, checkUpsolveFinished_StateMachine")
  if cliqst == :upsolved
      frsyms = getCliqFrontalVarIds(csmc.cliq)
    infocsm(csmc, "9, checkUpsolveFinished_StateMachine -- going for transferUpdateSubGraph! on $frsyms")
    # TODO what about down solve??
    transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frsyms, csmc.logger, updatePPE=false)

    # remove any solvable upward cached data -- TODO will have to be changed for long down partial chains
    # assuming maximally complte up solved cliq at this point
    # lockUpStatus!(csmc.cliq, csmc.cliq.index, true, csmc.logger, true, "9.finishCliqSolveCheck")
    sdims = Dict{Symbol,Float64}()
    for varid in getCliqAllVarIds(csmc.cliq)
      sdims[varid] = 0.0
    end
    updateCliqSolvableDims!(csmc.cliq, sdims, csmc.logger)
    # unlockUpStatus!(csmc.cliq)

    # go to 10
    return canCliqDownSolve_StateMachine # IncrementalInference.exitStateMachine
  elseif cliqst == :initialized
    setCliqDrawColor(csmc.cliq, "sienna")

    # go to 7
    return determineCliqNeedDownMsg_StateMachine
  else
    infocsm(csmc, "9, checkUpsolveFinished_StateMachine -- init not complete and should wait on init down message.")
    # setCliqDrawColor(csmc.cliq, "coral")
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

Do up initialization calculations,  preparation for solving Chapman-Kolmogorov
transit integral in upward direction.

Notes
- State machine function nr. 8f
- Includes initialization routines.
- Adds `:LIKELIHOODMESSAGE` factors but does not remove.
- gets msg likelihoods from cliqSubFg, see #760

DevNotes
- TODO: Make multi-core
"""
function prepInitUp_StateMachine(csmc::CliqStateMachineContainer)
  setCliqDrawColor(csmc.cliq, "green")

  # check if init is required and possible
  infocsm(csmc, "8f, prepInitUp_StateMachine -- going for doCliqAutoInitUpPart1!.")
  # get incoming clique up messages
  upmsgs = getMsgsUpInitChildren(csmc, skip=[csmc.cliq.index;])
  # Filter for usable messages
  ## FIXME joint decomposition as differential likelihoods conversion must still be done for init
  dellist = []
  for (chid, lm) in upmsgs
    if !(lm.status in [:initialized;:upsolved;:marginalized;:downsolved;:uprecycled])
      push!(dellist, chid)
    end
  end
  dellist .|> x->delete!(upmsgs, x)

  # remove all lingering upmessage likelihoods
  oldTags = lsf(csmc.cliqSubFg, tags=[:LIKELIHOODMESSAGE;])
  0 < length(oldTags) ? @warn("stale LIKELIHOODMESSAGE tags present in prepInitUp_StateMachine") : nothing
  oldFcts = oldTags .|> x->getFactor(csmc.cliqSubFg, x)
  # add incoming up messages as priors to subfg
  infocsm(csmc, "8f, prepInitUp_StateMachine -- adding up message factors")
  deleteMsgFactors!(csmc.cliqSubFg, oldFcts)
  # interally adds :LIKELIHOODMESSAGE, :UPWARD_DIFFERENTIAL, :UPWARD_COMMON to each of the factors
  msgfcts = addMsgFactors!(csmc.cliqSubFg, upmsgs, UpwardPass)

  # store the cliqSubFg for later debugging
  opts = getSolverParams(csmc.dfg)
  if opts.dbg
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforeupsolve"))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforeupsolve.pdf"))
  end

  # go to 8m
  return tryInitCliq_StateMachine
end


"""
    $SIGNATURES

Calculate the full upward Chapman-Kolmogorov transit integral solution approximation (i.e. upsolve).

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

  # Send upward message, NOTE consolidation WIP #459
  infocsm(csmc, "8g, doCliqUpSolveInitialized_StateMachine -- setting up messages with status = :upsolved")
  prepPutCliqueStatusMsgUp!(csmc, :upsolved)

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
  return checkUpsolveFinished_StateMachine
end



"""
$SIGNATURES

WIP to resolve 459 dwnMsg consolidation.  This is partly doing some kind of downsolve but seems out of place.

Notes
- State machine function nr. 10a

DevNotes
- FIXME, resolve/consolidate whatever is going on here
"""
function wipRedirect459Dwn_StateMachine(csmc::CliqStateMachineContainer)
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
  infocsm(csmc, "10a, wipRedirect459Dwn_StateMachine, parent status=$prntst.")
  if prntst != :downsolved
    infocsm(csmc, "10a, wipRedirect459Dwn_StateMachine, going around again.")
    return canCliqDownSolve_StateMachine
  end

  infocsm(csmc, "10a, wipRedirect459Dwn_StateMachine, going for down solve.")
  # go to 11
  return doCliqDownSolve_StateMachine
end


"""
    $SIGNATURES

Need description for this???

Notes
- State machine function nr. 8c

DevNotes
- Must consolidate as part of #459
"""
function waitChangeOnParentCondition_StateMachine(csmc::CliqStateMachineContainer)
  #
  # setCliqDrawColor(csmc.cliq, "coral")

  prnt = getParent(csmc.tree, csmc.cliq)
  if 0 < length(prnt)
    infocsm(csmc, "8c, waitChangeOnParentCondition_StateMachine, wait on parent=$(prnt[1].index) for condition notify.")
    # @sync begin
    #   @async begin
    #     sleep(1)
    #     notify(getSolveCondition(prnt[1]))
    #   end
      # wait but don't clear what is in the Condition (guess)
      wait(getSolveCondition(prnt[1]))
    # end
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

Do down solve calculations, loosely translates to solving Chapman-Kolmogorov
transit integral in downward direction.

Notes
- State machine function nr. 8e.ii.
  - Follows routines in 8c.
- Pretty major repeat of functionality, FIXME
- TODO: Make multi-core

DevNotes
- TODO Lots of cleanup required, especially from calling function.
- TODO move directly into a CSM state function
- TODO figure out the diference, redirection, ie relation to 8m??


CONSOLIDATED OLDER FUNCTION DOCS:
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
  # # TODO: is status of notify required here? either up or down msg??
  setCliqueStatus!(csmc.cliq, cliqst)

  # go to 8l
  return rmMsgLikelihoodsAfterDwn_StateMachine
end



"""
    $SIGNATURES

Check if the clique is fully initialized.

Notes
- State machine function nr. 8m

DevNotes
- TODO figure out the relation or possible consolidation with 8e.ii
"""
function tryInitCliq_StateMachine(csmc::CliqStateMachineContainer)
  # attempt initialize if necessary
  someInit = false
  if !areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq)
    # structure for all up message densities computed during this initialization procedure.
    varorder = getCliqVarInitOrderUp(csmc.tree, csmc.cliq)
    someInit = cycleInitByVarOrder!(csmc.cliqSubFg, varorder, logger=csmc.logger)
    # is clique fully upsolved or only partially?
    # print out the partial init status of all vars in clique
    printCliqInitPartialInfo(csmc.cliqSubFg, csmc.cliq, csmc.logger)
    infocsm(csmc, "8m, tryInitCliq_StateMachine -- someInit=$someInit, varorder=$varorder")
  end

  chldneed = doAnyChildrenNeedDwnMsg(getChildren(csmc.tree, csmc.cliq))
  allvarinit = areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq)
  infocsm(csmc, "8m, tryInitCliq_StateMachine -- someInit=$someInit, chldneed=$chldneed, allvarinit=$allvarinit")

  # redirect if any children needdownmsg
  if someInit || chldneed
    # prep down init message
    prepPutCliqueStatusMsgDwn!(csmc, :initialized)
    # # go to 7b
    # return slowIfChildrenNotUpSolved_StateMachine

    # go to 7e
    return slowWhileInit_StateMachine

    # (short cut) check again if all cliq vars have been initialized so that full inference can occur on clique
  elseif allvarinit
    infocsm(csmc, "8m, tryInitCliq_StateMachine -- all initialized")
    # TODO should this set status=:initialized? and notify???
    setCliqueStatus!(csmc.cliq, :initialized)

    # go to 8g.
    return doCliqUpSolveInitialized_StateMachine
  end

  infocsm(csmc, "8m, tryInitCliq_StateMachine -- not able to init all")
  # TODO Simplify this
  status = getCliqueStatus(csmc.cliq)

  if !(status == :initialized || length(getParent(csmc.tree, csmc.cliq)) == 0)
    # notify of results (big part of #459 consolidation effort)
    setCliqDrawColor(csmc.cliq, "orchid")
    prepPutCliqueStatusMsgUp!(csmc, :needdownmsg)
  end

  # go to 8h
  return rmUpLikeliSaveSubFg_StateMachine
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

Function to iterate through while initializing various child cliques that start off `needdownmsg`.

Notes
- State machine function nr. 7e
"""
function slowWhileInit_StateMachine(csmc::CliqStateMachineContainer)

  if doAnyChildrenNeedDwnMsg(csmc.tree, csmc.cliq)
    infocsm(csmc, "7e, slowWhileInit_StateMachine, must wait for new child messages.")

    # wait for THIS clique to be notified (PUSH NOTIFICATION FROM CHILDREN at `prepPutCliqueStatusMsgUp!`)
    wait(getSolveCondition(csmc.cliq))
  end

  # go to 8f
  return prepInitUp_StateMachine
end

"""
    $SIGNATURES

Delay loop if waiting on upsolves to complete.

Notes
- State machine 7b
- Differs from 4e in that here children must be "upsolved" or equivalent to continue.
  - Also waits on condition, so less succeptible to stale messages from children

DevNotes
- CONSIDER "recursion", to come back to this function to make sure that child clique updates are in fact upsolved.
  - other steps in CSM require each CSM to progress until cascading init problem is solved
"""
function slowIfChildrenNotUpSolved_StateMachine(csmc::CliqStateMachineContainer)

  # childs = getChildren(csmc.tree, csmc.cliq)
  # len = length(childs)
  for chld in getChildren(csmc.tree, csmc.cliq)
    chst = getCliqueStatus(chld)
    if !(chst in [:upsolved;:uprecycled;:marginalized;])
      infocsm(csmc, "7b, slowIfChildrenNotUpSolved_StateMachine, wait $(chst), cliq=$(chld.index), ch_lbl=$(getCliqFrontalVarIds(chld)[1]).")


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


"""
    $SIGNATURES

Block until all children have a csm status, using `fetch---Condition` model.

Notes
- State machine function nr.4e
- Blocking call if there are no status messages available
  - Will not block on stale message
"""
function blockUntilChildrenHaveStatus_StateMachine(csmc::CliqStateMachineContainer)
  #must happen before if :null

  notsolved = true

  while notsolved
    notsolved = false
    infocsm(csmc, "4e, blockUntilChildrenHaveStatus_StateMachine, get new status msgs.")

    stdict = fetchChildrenStatusUp(csmc.tree, csmc.cliq, csmc.logger)
    for (cid, st) in stdict
      infocsm(csmc, "4e, blockUntilChildrenHaveStatus_StateMachine, maybe wait cliq=$(cid), child status=$(st).")

      if st in [:null;]
        infocsm(csmc, "4e, blockUntilChildrenHaveStatus_StateMachine, waiting cliq=$(cid), child status=$(st).")
        wait(getSolveCondition(csmc.tree.cliques[cid]))
        notsolved = true
        break
      end
    end
  end
  
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
- FIXME, this guy is never called during WIP of 459 dwnMsg consolidation
- FIXME, something wrong with CSM sequencing, https://github.com/JuliaRobotics/IncrementalInference.jl/issues/602#issuecomment-682114232
- This might be replaced with 4-stroke tree-init.
"""
function dwnInitSiblingWaitOrder_StateMachine(csmc::CliqStateMachineContainer)

  prnt_ = getParent(csmc.tree, csmc.cliq)
  if 0 == length(prnt_)
    # go to 7
    return determineCliqNeedDownMsg_StateMachine
  end
  prnt = prnt_[1]

  opt = getSolverParams(csmc.cliqSubFg) # csmc.dfg

  # now get the newly computed message from the appropriate container
  # make sure this is a pull model #674 (pull msg from source/prnt)
  # FIXME must be consolidated as part of #459
  # dwinmsgs = getfetchCliqueInitMsgDown(getCliqueData(prnt), from=:dwnInitSiblingWaitOrder_StateMachine)
  dwinmsgs = fetchDwnMsgConsolidated(prnt)
  dwnkeys_ = collect(keys(dwinmsgs.belief))
  infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, dwinmsgs keys=$(dwnkeys_)")

  # add downward belief prop msgs
  msgfcts = addMsgFactors!(csmc.cliqSubFg, dwinmsgs, DownwardPass)
  ## update solveable dims with newly avaiable init information
  # determine if more info is needed for partial
  sdims = getCliqVariableMoreInitDims(csmc.cliqSubFg, csmc.cliq)
  infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, sdims=$(sdims)")
  updateCliqSolvableDims!(csmc.cliq, sdims, csmc.logger)

  if opt.dbg
    saveDFG(joinLogPath(csmc.cliqSubFg, "logs", "cliq$(csmc.cliq.index)", "fg_DWNCMN_8j"), csmc.cliqSubFg)
    drawGraph(csmc.cliqSubFg, show=false, filepath=(joinLogPath(csmc.cliqSubFg, "logs", "cliq$(csmc.cliq.index)", "fg_DWNCMN_8j.pdf")))
  end


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
    if dwinmsgs.status == :initialized && getCliqueStatus(csmc.cliq) == :needdowninit
      # go to 7e
      return slowWhileInit_StateMachine
    end

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
      return waitChangeOnParentCondition_StateMachine
    # # go to 8a
    # return collectDwnInitMsgFromParent_StateMachine
  # HALF DUPLICATED IN STEP 4
  elseif cliqst == :marginalized
    # go to 1
    return isCliqUpSolved_StateMachine
  end

  # go to 8b
  return attemptCliqInitUp_StateMachine
end


# """
#     $SIGNATURES

# Blocking case when all siblings and parent :needdownmsg.

# Notes
# - State machine function nr. 6c

# DevNotes
# - FIXME understand if this should be consolidated with 4b. `trafficRedirectConsolidate459_StateMachine`?
# - FIXME understand if this should be consolidated with 7c. `towardUpOrDwnSolve_StateMachine`
# """
# function doesParentNeedDwn_StateMachine(csmc::CliqStateMachineContainer)

#   infocsm(csmc, "6c, check/block sibl&prnt :needdownmsg")
#   prnt = getParent(csmc.tree, csmc.cliq)
#   if 0 == length(prnt) || getCliqueStatus(prnt[1]) != :needdownmsg
#     infocsm(csmc, "6c, prnt $(getCliqueStatus(prnt[1]))")
#     # go to 7
#     return determineCliqNeedDownMsg_StateMachine
#   # elseif
#   end

#   # go to 6d
#   return doAllSiblingsNeedDwn_StateMachine
# end


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


  # Some traffic direction
  if cliqst == :null
    # go to 4d
    return maybeNeedDwnMsg_StateMachine
  end

  # go to 6c
  # return doesParentNeedDwn_StateMachine
  prnt = getParent(csmc.tree, csmc.cliq)
  if 0 < length(prnt) && cliqst == :needdownmsg
    if getCliqueStatus(prnt[1]) == :needdownmsg
      # go to 8c
      return waitChangeOnParentCondition_StateMachine
    end
    # go to 6d
    return doAllSiblingsNeedDwn_StateMachine
  end

  # go to 7
  return determineCliqNeedDownMsg_StateMachine
end




## ============================================================================================
# START of dwnmsg consolidation bonanza
## ============================================================================================
#
# blockCliqSiblingsParentChildrenNeedDown_  # Blocking case when all siblings and parent :needdownmsg.
# doAllSiblingsNeedDwn_StateMachine   # Trying to figure out when to block on siblings for cascade down init.
# maybeNeedDwnMsg_StateMachine         # If all children (then also escalate to) :needdownmsgs and block until sibling status.
# determineCliqNeedDownMsg_StateMachine     # Try decide whether this `csmc.cliq` needs a downward initialization message.
# doAnyChildrenNeedDwn_StateMachine         # Determine if any one of the children :needdownmsg.
# downInitRequirement_StateMachine          # Place fake up msg and notify down init status if any children :needdownmsg




"""
$SIGNATURES

If all children (then also escalate to) :needdownmsgs and block until sibling status.

Notes
- State machine function nr.4d

DevNotes
- TODO consolidate with 6d?????
- TODO Any overlap with nr.4c??
"""
function maybeNeedDwnMsg_StateMachine(csmc::CliqStateMachineContainer)

  # fetch (should not block)
  stdict = fetchChildrenStatusUp(csmc.tree, csmc.cliq, csmc.logger)
  chstatus = collect(values(stdict))
  infocsm(csmc,"fetched all, keys=$(keys(stdict)), values=$(chstatus).")
  len = length(chstatus)

  # if all children needdownmsg
  if 0 < len && sum(chstatus .== :needdownmsg) == len
    # Can this cliq init with local information?
    # get initial estimate of solvable dims in this cliq
    sdims = getCliqVariableMoreInitDims(csmc.cliqSubFg, csmc.cliq)
    infocsm(csmc, "4d, maybeNeedDwnMsg_StateMachine, sdims=$(sdims)")
    updateCliqSolvableDims!(csmc.cliq, sdims, csmc.logger)
    if 0 < sum(collect(values(sdims)))
      # try initialize
      # go to 8f
      return prepInitUp_StateMachine
    end

    # TODO maybe can happen where some children need more information?
    infocsm(csmc, "4d, maybeNeedDwnMsg_StateMachine, escalating to :needdownmsg since all children :needdownmsg")
    # NOTE, trying consolidation with prepPutUp for #459 effort
    setCliqDrawColor(csmc.cliq, "orchid1")
    prepPutCliqueStatusMsgUp!(csmc, :needdownmsg)

    # debuggin #459 transition
    infocsm(csmc, "4d, maybeNeedDwnMsg_StateMachine -- finishing before going to  blockSiblingStatus_StateMachine")

    # go to 5
    return blockSiblingStatus_StateMachine
  end

  if doAnyChildrenNeedDwnMsg(csmc.tree, csmc.cliq)
    # go to 7e (#754)
    return slowWhileInit_StateMachine
  end

  # go to 7
  return determineCliqNeedDownMsg_StateMachine
end


"""
    $SIGNATURES

Try decide whether this `csmc.cliq` needs a downward initialization message.

Notes
- State machine function nr. 7

DevNotes
- Consolidate with 4b?
- found by trail and error, TODO review and consolidate with rest of CSM after major #459 and PCSM consolidation work is done.
"""
function determineCliqNeedDownMsg_StateMachine(csmc::CliqStateMachineContainer)

  # fetch children status
  stdict = fetchChildrenStatusUp(csmc.tree, csmc.cliq, csmc.logger)
  infocsm(csmc,"fetched all, keys=$(keys(stdict)).")

  # hard assumption here on upsolve from leaves to root
  childst = collect(values(stdict))
  # are child cliques sufficiently solved
  resolveinit = (filter(x-> x in [:upsolved;:marginalized;:downsolved;:uprecycled], childst) |> length) == length(childst)
  chldupandinit = sum(childst .|> x-> (x in [:initialized;:upsolved;:marginalized;:downsolved;:uprecycled])) == length(childst)
  allneeddwn = (filter(x-> x == :needdownmsg, childst) |> length) == length(childst) && 0 < length(childst)
  chldneeddwn = :needdownmsg in childst
  chldnull = :null in childst

  cliqst = getCliqueStatus(csmc.cliq)

  infocsm(csmc, "7, determineCliqNeedDownMsg_StateMachine, childst=$childst")
  infocsm(csmc, "7, determineCliqNeedDownMsg_StateMachine, cliqst=$cliqst, resolveinit=$resolveinit, allneeddwn=$allneeddwn, chldneeddwn=$chldneeddwn, chldupandinit=$chldupandinit, chldnull=$chldnull")


  # merged in from 4c here into 7, part of dwnMsg #459
  if cliqst == :needdownmsg
    if allneeddwn || length(stdict) == 0
      infocsm(csmc, "7, determineCliqNeedDownMsg_StateMachine, at least some children :needdownmsg")
      # # go to 8j
      return dwnInitSiblingWaitOrder_StateMachine
    end
    # FIXME FIXME FIXME hex 2.16 & 3.21 should go for downinitfrom parent
    # perhaps should check if all siblings also needdownmsg???

    if chldneeddwn
      # go to 5
      return blockSiblingStatus_StateMachine
      # # go to 7b
      # return slowIfChildrenNotUpSolved_StateMachine
    end
  end

  # includes case of no children
  if resolveinit && !chldneeddwn
    # go to 8j (dwnMsg #459 WIP 9)
    # return dwnInitSiblingWaitOrder_StateMachine
    # go to 7c
    return towardUpOrDwnSolve_StateMachine
  end

  if chldnull || cliqst == :null && !chldupandinit
    # go to 4e
    return blockUntilChildrenHaveStatus_StateMachine
  elseif chldneeddwn || chldupandinit
    # go to 7b
    return slowIfChildrenNotUpSolved_StateMachine
  end

  # go to 6d
  return doAllSiblingsNeedDwn_StateMachine
end


"""
    $SIGNATURES

Detect block on clique's parent.

Notes
- State machine function nr. 5

DevNotes
- FIXME refactor this for when prnt==:null, cliq==:needdowninit
"""
function blockSiblingStatus_StateMachine(csmc::CliqStateMachineContainer)
    # infocsm(csmc, "5, blocking on parent until all sibling cliques have valid status")
    # setCliqDrawColor(csmc.cliq, "blueviolet")

    cliqst = getCliqueStatus(csmc.cliq)
    infocsm(csmc, "5, block on siblings")
    prnt = getParent(csmc.tree, csmc.cliq)
    if cliqst == :needdownmsg && 0 < length(prnt) && getCliqueStatus(prnt[1]) in [:null;:needdownmsg]
      # go to 8c
      return waitChangeOnParentCondition_StateMachine
    end

    # infocsm(csmc, "5, has parent clique=$(prnt[1].index)")
    # ret = fetchChildrenStatusUp(csmc.tree, prnt[1], csmc.logger)
    # infocsm(csmc,"prnt $(prnt[1].index), fetched all, keys=$(keys(ret)).")

  infocsm(csmc, "5, finishing")
  # go to 6c
  # return doesParentNeedDwn_StateMachine

  # go to 6d
  return doAllSiblingsNeedDwn_StateMachine
end


"""
$SIGNATURES

Trying to figure out when to block on siblings for cascade down init.

Notes
- State machine function nr.6d
- Assume there must be a parent.
- Part of #459 dwnMsg consolidation work.
- used for regulating long need down message chains.
- exit strategy is parent becomes status `:initialized`.

DevNotes
- Consolidation with work with similar likely required.
"""
function doAllSiblingsNeedDwn_StateMachine(csmc::CliqStateMachineContainer)

  # go to 6c
  # return doesParentNeedDwn_StateMachine
  infocsm(csmc, "7, check/block sibl&prnt :needdownmsg")
  prnt = getParent(csmc.tree, csmc.cliq)
  if 0 == length(prnt) # || getCliqueStatus(prnt[1]) != :needdownmsg
    # infocsm(csmc, "4d, prnt $(getCliqueStatus(prnt[1]))")
    # go to 7
    return determineCliqNeedDownMsg_StateMachine
  # elseif
  end


  prnt = getParent(csmc.tree, csmc.cliq)
  prnt_ = prnt[1]

  stdict = fetchChildrenStatusUp(csmc.tree, prnt_, csmc.logger)
  # remove this clique's data from dict
  delete!(stdict, csmc.cliq.index)
  # hard assumption here on upsolve from leaves to root
  siblst = values(collect(stdict))
  # are child cliques sufficiently solved
  allneeddwn = (filter(x-> x == :needdownmsg, siblst) |> length) == length(siblst) && 0 < length(siblst)


  # FIXME, understand why is there another status event from parent msg here... How to consolidate this with CSM 8a
  if allneeddwn
    # # go to 6e
    # return slowOnPrntAsChildrNeedDwn_StateMachine
    # go to 8j
    return dwnInitSiblingWaitOrder_StateMachine
  end

  # go to 7
  return determineCliqNeedDownMsg_StateMachine
end


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
  # for ch in getChildren(csmc.tree, csmc.cliq)
  for (clid, chst) in fetchChildrenStatusUp(csmc.tree, csmc.cliq, csmc.logger)
    # if getCliqueStatus(ch) == :needdownmsg
    if chst == :needdowninit
      someChildrenNeedDwn = true
      break
    end
  end

  if someChildrenNeedDwn
    # send a down init message
    prepPutCliqueStatusMsgDwn!(csmc)
    # go to 8k
    return sendCurrentUpMsg_StateMachine
  end

  # go to 8b
  return attemptCliqInitUp_StateMachine
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
    return prepInitUp_StateMachine
  end

  if !notChildNeedDwn
    # go to 8m
    return tryInitCliq_StateMachine
  end

  # go to 9
  return checkUpsolveFinished_StateMachine
end


## ============================================================================================
# End of dwnmsg consolidation bonanza
## ============================================================================================




#
