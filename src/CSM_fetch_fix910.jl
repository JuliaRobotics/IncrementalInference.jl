# Downward initialization through a cascading model
# these functions need to be consolidated refactored and deprecated accordingly
# see #910


## ============================================================================================
## Down init sibling priority order
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

  _dbgCSMSaveSubFG(csmc, "fg_DWNCMN_8j")

  ## FIXME, if this new, and all sibling clique's solvableDim are 0, then go back to waitChangeOnParentCondition_StateMachine

  # go to 8o.i
  return testDirectDwnInit_StateMachine
end

"""
    $SIGNATURES

Can this clique initialize directly from available down message info?

Notes
- State machine function nr. 8o.i
- Assume must have parent since this only occurs after 8j.
"""
function testDirectDwnInit_StateMachine(csmc::CliqStateMachineContainer)

  prnt = getParent(csmc.tree, csmc.cliq)[1]
  dwinmsgs = fetchDwnMsgConsolidated(prnt)
  dwnkeys_ = collect(keys(dwinmsgs.belief))
  # NOTE, only use separators, not all parent variables (DF ???)
  # dwnkeys_ = lsf(csmc.cliqSubFg, tags=[:DOWNWARD_COMMON;]) .|> x->ls(csmc.cliqSubFg, x)[1]
  # @assert length(intersect(dwnkeys, dwnkeys_)) == length(dwnkeys) "split dwnkeys_ is not the same, $dwnkeys, and $dwnkeys_, separators: $(getCliqSeparatorVarIds(csmc.cliq))"

  # priorize solve order for mustinitdown with lowest dependency first
  # follow example from issue #344
  # mustwait = false
  if length(intersect(dwnkeys_, getCliqSeparatorVarIds(csmc.cliq))) == 0
    # can directly use DOWNWARD_COMMON
    infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, no can do, must wait for siblings to update parent first.")
    # mustwait = true

    # go to 8o.iv
    return sibsDwnPriorityInit_StateMachine
  end

  # go to 8o.ii
  return testDelayOrderDwnInit_StateMachine
end

"""
    $SIGNATURES


Notes
- State machine function nr. 8o.ii

"""
function testDelayOrderDwnInit_StateMachine(csmc::CliqStateMachineContainer)

  prnt = getParent(csmc.tree, csmc.cliq)[1]
  dwinmsgs = fetchDwnMsgConsolidated(prnt)
  dwnkeys_ = collect(keys(dwinmsgs.belief))
  
  if getSiblingsDelayOrder(csmc.tree, csmc.cliq, dwnkeys_, logger=csmc.logger)
    # check sibling status and decide where next based on available dwn msgs
    infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, prioritize")
    # mustwait = true

    # go to 8o.iv
    return sibsDwnPriorityInit_StateMachine
  end

  # go to 8o.iii
  return testPartialNeedsDwnInit_StateMachine
end

"""
    $SIGNATURES


Notes
- State machine function nr. 8o.iii
"""
function testPartialNeedsDwnInit_StateMachine(csmc::CliqStateMachineContainer)
  prnt = getParent(csmc.tree, csmc.cliq)[1]
  dwinmsgs = fetchDwnMsgConsolidated(prnt)
  
  if getCliqSiblingsPartialNeeds(csmc.tree, csmc.cliq, dwinmsgs, logger=csmc.logger)
    # if both requires and more down msg info is likely
    infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, partialneedsmore")
    # mustwait = true

    # go to 8o.iv
    return sibsDwnPriorityInit_StateMachine
  end

  # go to 8e.ii.
  return tryDwnInitCliq_StateMachine
end

"""
    $SIGNATURES


Notes
- State machine function nr. 8o.iv
"""
function sibsDwnPriorityInit_StateMachine(csmc::CliqStateMachineContainer)

  prnt = getParent(csmc.tree, csmc.cliq)[1]
  dwinmsgs = fetchDwnMsgConsolidated(prnt)
  dwnkeys_ = collect(keys(dwinmsgs.belief))

  solord = getCliqSiblingsPriorityInitOrder( csmc.tree, prnt, dwinmsgs, csmc.logger )
  noOneElse = areSiblingsRemaingNeedDownOnly(csmc.tree, csmc.cliq)
  infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, $(prnt.index), $noOneElse, solord = $solord")


  if csmc.cliq.index != solord[1]
    # TODO, is this needed? fails hex if included with correction :needdowninit --> :needdownmsg
    # if dwinmsgs.status == :initialized && getCliqueStatus(csmc.cliq) == :needdownmsg  # used to be :needdowninit
    #   # go to 7e
    #   return slowWhileInit_StateMachine
    # end

    infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, must wait on change.")
    # remove all message factors
    # remove msg factors previously added
    fctstorm = deleteMsgFactors!(csmc.cliqSubFg, [:DOWNWARD_COMMON])
    infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, removing factors $fctstorm")

    # go to 8c
    return waitChangeOnParentCondition_StateMachine
  end

  # go to 8e.ii.
  return tryDwnInitCliq_StateMachine
end








