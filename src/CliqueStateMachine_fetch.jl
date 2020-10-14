# clique state machine for tree based initialization and inference

# newer exports
# export towardUpOrDwnSolve_StateMachine, maybeNeedDwnMsg_StateMachine
# export prepInitUp_StateMachine, doCliqUpSolveInitialized_StateMachine
# export rmUpLikeliSaveSubFg_StateMachine
# export blockCliqSiblingsParentChildrenNeedDown_StateMachine

## ============================================================================================
## Only fetch related functions after this, all to be deprecated
## ============================================================================================


# XXX only fetch
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


# XXX only fetch
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



# XXX only fetch
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





# XXX only fetch
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


# XXX only fetch
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




# XXX only fetch
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


# XXX only fetch
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


# XXX only fetch
"""
    $SIGNATURES

Detect block on clique's parent.

Notes
- State machine function nr. 5

DevNotes
- FIXME refactor this for when prnt==:null, cliq==:needdownmsg/init
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


# XXX only fetch
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


# XXX only fetch
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
    if chst == :needdownmsg # NOTE was :needdowninit
      someChildrenNeedDwn = true
      break
    end
  end

  if someChildrenNeedDwn
    # send a down init message
    # TODO down message is sent more than once?
    # here and in sendCurrentUpMsg_StateMachine
    prepPutCliqueStatusMsgDwn!(csmc)
    # go to 8k
    return sendCurrentUpMsg_StateMachine
  end

  # go to 8b
  return attemptCliqInitUp_StateMachine
end


# XXX only fetch
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
    return tryUpInitCliq_StateMachine
  end

  # go to 9
  return checkUpsolveFinished_StateMachine
end


## ============================================================================================
# End of dwnmsg consolidation bonanza
## ============================================================================================




#
