## ==================================================================================================================='
# Should be made Common CSM functions
## ==================================================================================================================='

"""
    $SIGNATURES

Notify possible parent if clique is upsolved and exit the state machine.

Notes
- State machine function nr.0
- can recycle if two checks:
  - previous clique was identically downsolved
  - all children are also :uprecycled
"""
function canCliqMargRecycle_StateMachine(csmc::CliqStateMachineContainer)
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
  return canCliqIncrRecycle_StateMachine
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
function canCliqIncrRecycle_StateMachine(csmc::CliqStateMachineContainer)
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

Build a sub factor graph for clique variables from the larger factor graph.

Notes
- State machine function nr.2
"""
function buildCliqSubgraph_StateMachine(csmc::CliqStateMachineContainer)
  # build a local subgraph for inference operations
  infocsm(csmc, "2, build subgraph syms=$(getCliqAllVarIds(csmc.cliq))")
  buildCliqSubgraph!(csmc.cliqSubFg, csmc.dfg, csmc.cliq)

  # if dfg, store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_build")

  # go to 4
  return canCliqMargSkipUpSolve_StateMachine
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
  # marginalized state is set in `canCliqMargRecycle_StateMachine`
  if cliqst == :marginalized
    # go to 10 -- Add case for IIF issue #474
    return canCliqDownSolve_StateMachine
  end

  # go to 4e
  return blockUntilChildrenHaveStatus_StateMachine
end


## ==================================================================================================================='
## Does this have a place
## ==================================================================================================================='

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
  _dbgCSMSaveSubFG(csmc, "fg_build_down")

  # go to 10
  return canCliqDownSolve_StateMachine
end


## ==================================================================================================================='
## Split and use
## ==================================================================================================================='

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
  _dbgCSMSaveSubFG(csmc, "fg_afterdownsolve")

  # transfer results to main factor graph
  frsyms = getCliqFrontalVarIds(csmc.cliq)
  infocsm(csmc, "11, finishingCliq -- going for transferUpdateSubGraph! on $frsyms")
  transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frsyms, csmc.logger, updatePPE=true)

  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- before prepPutCliqueStatusMsgDwn!")
  cliqst = prepPutCliqueStatusMsgDwn!(csmc, :downsolved)
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- just notified prepPutCliqueStatusMsgDwn!")

  # remove msg factors that were added to the subfg
  rmFcts = deleteMsgFactors!(csmc.cliqSubFg)
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- removing all up/dwn message factors, length=$(length(rmFcts))")

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
  putDwnMsgConsolidated!(csmc.cliq.data, dwnmsgs) # , from=:putMsgDwnThis!  putCliqueMsgDown!
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

  prepPutCliqueStatusMsgDwn!(csmc, :downsolved)
  # notifyCliqDownInitStatus!(csmc.cliq, :downsolved, logger=csmc.logger)

  # bye
  return IncrementalInference.exitStateMachine
end




## ==================================================================================================================='
## Old 
## ==================================================================================================================='

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
  upmsg = prepCliqueMsgUpConsolidated(csmc.cliqSubFg, csmc.cliq, getCliqueStatus(csmc.cliq), logger=csmc.logger)
  prepPutCliqueStatusMsgUp!(csmc, upmsg=upmsg)

  # also send a down message -- seem weird while doing #459 but okay
  cliqst = prepPutCliqueStatusMsgDwn!(csmc)

  # Legend: initialized but not solved yet (likely child cliques that depend on downward autoinit msgs),
  setCliqDrawColor(csmc.cliq, "sienna")

  infocsm(csmc, "8k, sendCurrentUpMsg_StateMachine -- near-end down init attempt, $cliqst.")

  # go to 8b
  return attemptCliqInitUp_StateMachine
end










## =================================================================================================================


"""
    $SIGNATURES

Special internal function to try return the clique data if succesfully identified in `othertree::AbstractBayesTree`,
based on contents of `seeksSimilar::BayesTreeNodeData`.

Notes
- Used to identify and skip similar cliques (i.e. recycle computations)
"""
function attemptTreeSimilarClique(othertree::AbstractBayesTree, 
                                  seeksSimilar::BayesTreeNodeData  )
  #
  # inner convenience function for returning empty clique
  function EMPTYCLIQ()
    clq = TreeClique(-1,"null")
    setLabel!(clq, "")
    setCliqueData!(clq, BayesTreeNodeData())
    return clq
  end

  # does the other clique even exist?
  seekFrontals = getCliqFrontalVarIds(seeksSimilar)
  if !hasClique(othertree, seekFrontals[1])
    return EMPTYCLIQ()
  end

  # do the cliques share the same frontals?
  otherCliq = getClique(othertree, seekFrontals[1])
  otherFrontals = getCliqFrontalVarIds(otherCliq)
  commonFrontals = intersect(seekFrontals, otherFrontals)
  if length(commonFrontals) != length(seekFrontals) || length(commonFrontals) != length(otherFrontals)
    return EMPTYCLIQ()
  end

  # do the cliques share the same separator variables?
  seekSeparator = getCliqSeparatorVarIds(seeksSimilar)
  otherSeparator = getCliqSeparatorVarIds(otherCliq)
  commonSep = intersect(seekSeparator, otherSeparator)
  if length(commonSep) != length(seekSeparator) || length(commonSep) != length(otherSeparator)
    return EMPTYCLIQ()
  end

  # do the cliques use the same factors (potentials)
  seekPotentials = getCliqFactorIds(seeksSimilar)
  otherFactors = getCliqFactorIds(otherCliq)
  commonFactors = intersect(seekPotentials, otherFactors)
  if length(commonFactors) != length(seekPotentials) || length(commonFactors) != length(otherFactors)
    return EMPTYCLIQ()
  end

  # lets assume they are the same
  return otherCliq::TreeClique
end



## ==================================================================================================================='
## Can be consolidated/used (mostly used already as copies in X)
## ==================================================================================================================='


"""
    $SIGNATURES

Do cliq downward inference

Notes:
- State machine function nr. 11
"""
function doCliqDownSolve_StateMachine(csmc::CliqStateMachineContainer)
  infocsm(csmc, "11, doCliqDownSolve_StateMachine")
  setCliqDrawColor(csmc.cliq, "red")
  # get down msg from parent (assuming root clique CSM wont make it here)
  # this looks like a pull model #674
  prnt = getParent(csmc.tree, csmc.cliq)
  dwnmsgs = fetchDwnMsgConsolidated(prnt[1])
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- dwnmsgs=$(collect(keys(dwnmsgs.belief)))")

  __doCliqDownSolve!(csmc, dwnmsgs)

  # compute new down messages
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- going to set new down msgs.")
  newDwnMsgs = getSetDownMessagesComplete!(csmc.cliqSubFg, csmc.cliq, dwnmsgs, csmc.logger)

  prepPutCliqueStatusMsgDwn!(csmc, :downsolved)

  # update clique subgraph with new status
  setCliqDrawColor(csmc.cliq, "lightblue")

  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- finished with downGibbsCliqueDensity, now update csmc")

  # go to 11b.
  return cleanupAfterDownSolve_StateMachine
end

# XXX only use skip down part
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

  # go to 8c
  return waitChangeOnParentCondition_StateMachine
  # # go to 10a
  # return wipRedirect459Dwn_StateMachine
end

# XXX does not look like it has a place
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
    # setCliqDrawColor(csmc.cliq, "sienna")

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

# XXX does not look like it has a place
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
  oldTags = deleteMsgFactors!(csmc.cliqSubFg)
  0 < length(oldTags) ? @warn("stale LIKELIHOODMESSAGE tags present in prepInitUp_StateMachine") : nothing
  
  # add incoming up messages as priors to subfg
  infocsm(csmc, "8f, prepInitUp_StateMachine -- adding up message factors")
  # interally adds :LIKELIHOODMESSAGE, :UPWARD_DIFFERENTIAL, :UPWARD_COMMON to each of the factors
  msgfcts = addMsgFactors!(csmc.cliqSubFg, upmsgs, UpwardPass)

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_beforeupsolve")

  # go to 8m
  return tryUpInitCliq_StateMachine
end

#XXX uses same functions, can be split at message
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

  __doCliqUpSolveInitialized!(csmc)

  # Send upward message, NOTE consolidation WIP #459
  infocsm(csmc, "8g, doCliqUpSolveInitialized_StateMachine -- setting up messages with status = :upsolved")
  prepPutCliqueStatusMsgUp!(csmc, :upsolved)

  # go to 8h
  return rmUpLikeliSaveSubFg_StateMachine
end


# XXX can be changed to utility function
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
  opts = getSolverParams(csmc.dfg)

  # remove msg factors that were added to the subfg
  tags__ = opts.useMsgLikelihoods ? [:UPWARD_COMMON;] : [:LIKELIHOODMESSAGE;]
  msgfcts = deleteMsgFactors!(csmc.cliqSubFg, tags__)
  infocsm(csmc, "8g, doCliqUpsSolveInit.! -- status = $(status), removing $(tags__) factors, length=$(length(msgfcts))")

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_afterupsolve")

  # go to 9
  return checkUpsolveFinished_StateMachine
end


# XXX does not look like it has a place, replaced by take!
"""
    $SIGNATURES

When this clique needs information from parent to continue but parent is still busy.  
Reasons are either downsolve or need down init message information (which are similar).

Notes
- State machine function nr. 8c
- bad idea to injectDelayBefore this function, because it will delay waiting on the parent past the event.
"""
function waitChangeOnParentCondition_StateMachine(csmc::CliqStateMachineContainer)
  #
  # setCliqDrawColor(csmc.cliq, "coral")

  prnt = getParent(csmc.tree, csmc.cliq)
  if 0 < length(prnt)
    infocsm(csmc, "8c, waitChangeOnParentCondition_StateMachine, wait on parent=$(prnt[1].index) for condition notify.")
    
    prntst = fetchDwnMsgConsolidated(prnt[1]).status
    if prntst != :downsolved
      wait(getSolveCondition(prnt[1]))
    end

    # wait(getSolveCondition(prnt[1]))
  else
    infocsm(csmc, "8c, waitChangeOnParentCondition_StateMachine, cannot wait on parent for condition notify.")
    @warn "no parent!"
  end

  # Listing most likely status values that might lead to 4b (TODO needs to be validated)
  if getCliqueStatus(csmc.cliq) in [:needdownmsg; :initialized; :null]
    # go to 4b
    return trafficRedirectConsolidate459_StateMachine
  end

  ## consolidation from CSM 10a
  # yes, continue with downsolve
  infocsm(csmc, "10a, wipRedirect459Dwn_StateMachine, parent status=$prntst.")
  if prntst != :downsolved
    infocsm(csmc, "10a, wipRedirect459Dwn_StateMachine, going around again.")
    return canCliqDownSolve_StateMachine
  end

  infocsm(csmc, "10a, wipRedirect459Dwn_StateMachine, going for down solve.")
  # go to 11
  return doCliqDownSolve_StateMachine
end

# XXX functions replaced in take! structure
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
function tryDwnInitCliq_StateMachine(csmc::CliqStateMachineContainer)
  setCliqDrawColor(csmc.cliq, "green")
  opt = getSolverParams(csmc.cliqSubFg)
  dwnkeys_ = lsf(csmc.cliqSubFg, tags=[:DOWNWARD_COMMON;]) .|> x->ls(csmc.cliqSubFg, x)[1]

  ## TODO deal with partial inits only, either delay or continue at end...
  # find intersect between downinitmsgs and local clique variables
  # if only partials available, then
  infocsm(csmc, "8e.ii, tryDwnInitCliq_StateMachine, do cliq init down dwinmsgs=$(dwnkeys_)")
  # get down variable initialization order
  initorder = getCliqInitVarOrderDown(csmc.cliqSubFg, csmc.cliq, dwnkeys_)
  infocsm(csmc, "8e.ii, tryDwnInitCliq_StateMachine,  initorder=$(initorder)")

  # store the cliqSubFg for later debugging
  if opt.dbg
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opt.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforedowninit"))
  end

  # cycle through vars and attempt init
  infocsm(csmc, "8e.ii, tryDwnInitCliq_StateMachine, cycle through vars and attempt init")
  # cliqst = :needdownmsg
  if cycleInitByVarOrder!(csmc.cliqSubFg, initorder, logger=csmc.logger)
    # cliqst = :initialized
    # TODO: transfer values changed in the cliques should be transfered to the tree in proc 1 here.
    # # TODO: is status of notify required here? either up or down msg??
    setCliqDrawColor(csmc.cliq, "sienna")
    setCliqueStatus!(csmc.cliq, :initialized)
  end
  
  # go to 8l
  return rmMsgLikelihoodsAfterDwn_StateMachine
end


# XXX functions replaced in take! structure
"""
    $SIGNATURES

Check if the clique is fully initialized.

Notes
- State machine function nr. 8m

DevNotes
- TODO figure out the relation or possible consolidation with 8e.ii
"""
function tryUpInitCliq_StateMachine(csmc::CliqStateMachineContainer)
  # attempt initialize if necessary
  setCliqDrawColor(csmc.cliq, "green")
  someInit = false
  if !areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq)
    # structure for all up message densities computed during this initialization procedure.
    varorder = getCliqVarInitOrderUp(csmc.tree, csmc.cliq)
    someInit = cycleInitByVarOrder!(csmc.cliqSubFg, varorder, logger=csmc.logger)
    # is clique fully upsolved or only partially?
    # print out the partial init status of all vars in clique
    printCliqInitPartialInfo(csmc.cliqSubFg, csmc.cliq, csmc.logger)
    infocsm(csmc, "8m, tryUpInitCliq_StateMachine -- someInit=$someInit, varorder=$varorder")
  end

  chldneed = doAnyChildrenNeedDwnMsg(getChildren(csmc.tree, csmc.cliq))
  allvarinit = areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq)
  infocsm(csmc, "8m, tryUpInitCliq_StateMachine -- someInit=$someInit, chldneed=$chldneed, allvarinit=$allvarinit")

  upmessages = fetchMsgsUpChildrenDict(csmc)
  all_child_status = map(msg -> msg.status, values(upmessages))

  # redirect if any children needdownmsg
  if someInit || chldneed
    # Calculate and share the children sum solvableDim information for priority initialization
    totSolDims = Dict{Int, Float64}()
    for (clid, upmsg) in upmessages
      totSolDims[clid] = 0
      for (varsym, tbup) in upmsg.belief
        totSolDims[clid] += tbup.solvableDim
      end
    end
    infocsm(csmc, "8m, tryUpInitCliq_StateMachine -- totSolDims=$totSolDims")

    # prep and put down init message
    setCliqDrawColor(csmc.cliq, "sienna")
    prepPutCliqueStatusMsgDwn!(csmc, :initialized, childSolvDims=totSolDims)

    # go to 7e
    return slowWhileInit_StateMachine

    # (short cut) check again if all cliq vars have been initialized so that full inference can occur on clique
  # clique should be initialized and all children upsolved, uprecycled, or marginalized
  elseif allvarinit && all(in.(all_child_status, Ref([:upsolved; :uprecycled; :marginalized])))
    infocsm(csmc, "8m, tryUpInitCliq_StateMachine -- all initialized")
    setCliqDrawColor(csmc.cliq, "sienna")
    # don't send a message yet since the upsolve is about to occur too
    setCliqueStatus!(csmc.cliq, :initialized)

    # go to 8g.
    return doCliqUpSolveInitialized_StateMachine
  end

  infocsm(csmc, "8m, tryUpInitCliq_StateMachine -- not able to init all")
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

# XXX maybe change to utility function
"""
    $SIGNATURES

Remove any `:LIKELIHOODMESSAGE` from `cliqSubFg`.

Notes
- State machine function nr.8l
"""
function rmMsgLikelihoodsAfterDwn_StateMachine(csmc::CliqStateMachineContainer)
  ## TODO only remove :DOWNWARD_COMMON messages here
  #

  _dbgCSMSaveSubFG(csmc, "fg_afterdowninit")

  ## FIXME move this to separate state in CSM.
  # remove all message factors
  # remove msg factors previously added
  fctstorm = deleteMsgFactors!(csmc.cliqSubFg)
  infocsm(csmc, "8e.ii., tryDwnInitCliq_StateMachine, removing factors $fctstorm")

  # go to 8d
  return decideUpMsgOrInit_StateMachine
end
