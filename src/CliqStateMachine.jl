# clique state machine for tree based initialization and inference

"""
    $SIGNATURES

Specialized info logger print function to show clique state machine information
in a standardized form.
"""
function infocsm(csmc::CliqStateMachineContainer, str::A) where {A <: AbstractString}

  tm = string(Dates.now())
  tmt = split(tm, 'T')[end]

  lbl = csmc.cliq.attributes["label"]
  lbl1 = split(lbl,',')[1]
  cliqst = getCliqStatus(csmc.cliq)

  @info "$tmt | $(current_task()) cliq $(csmc.cliq.index), $lbl1, $(cliqst) -- "*str
  nothing
end


"""
    $SIGNATURES

Notes
- State machine function nr.9
"""
function finishCliqSolveCheck_StateMachine(csmc::CliqStateMachineContainer)
  csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing
  cliqst = getCliqStatus(csmc.cliq)
  infocsm(csmc, "9, finishingCliq")
  if cliqst == :upsolved
    infocsm(csmc, "9, going for transferUpdateSubGraph!")
    frsyms = getCliqFrontalVarIds(csmc.cliq)
    transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frsyms)
    return IncrementalInference.exitStateMachine
  elseif cliqst == :initialized
    setCliqDrawColor(csmc.cliq, "sienna")
    csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing
    # go to 7
    return determineCliqNeedDownMsg_StateMachine
  else
    infocsm(csmc, "9, init not complete and should wait on init down message.")
    setCliqDrawColor(csmc.cliq, "green")
    # TODO, potential problem with trying to downsolve
    # return isCliqNull_StateMachine # doesCliqNeeddownmsg_StateMachine
  end

  # go to 4
  return isCliqNull_StateMachine # whileCliqNotSolved_StateMachine
end

"""
    $SIGNATURES

Notes
- State machine function nr. 8c
"""
function waitChangeOnParentCondition_StateMachine(csmc::CliqStateMachineContainer)
  prnt = getParent(csmc.tree, csmc.cliq)
  if length(prnt) > 0
    infocsm(csmc, "waitChangeOnParentCondition_StateMachine, wait on parent for condition notify.")
    wait(getSolveCondition(prnt[1]))
  else
    infocsm(csmc, "waitChangeOnParentCondition_StateMachine, cannot wait on parent for condition notify.")
    @warn "no parent!"
  end

  # go to 4
  return isCliqNull_StateMachine
end


"""
    $SIGNATURES

Do up initialization calculations, loosely translates to solving Chapman-Kolmogorov
transit integral in upward direction.

Notes
- State machine function nr. 8b
- Includes initialization routines.
- TODO: Make multi-core
"""
function attemptCliqInitUp_StateMachine(csmc::CliqStateMachineContainer)


  cliqst = getCliqStatus(csmc.cliq)

  infocsm(csmc, "8b, doCliqAutoInitUp, !areCliqChildrenNeedDownMsg()=$(!areCliqChildrenNeedDownMsg(csmc.tree, csmc.cliq))" )
  if cliqst in [:initialized; :null; :needdownmsg] && !areCliqChildrenNeedDownMsg(csmc.tree, csmc.cliq)
    setCliqDrawColor(csmc.cliq, "red")
    csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing
    cliqst = doCliqAutoInitUp!(csmc.cliqSubFg, csmc.tree, csmc.cliq)
  end

  # go to 9
  return finishCliqSolveCheck_StateMachine
end


"""
    $SIGNATURES

Do down initialization calculations, loosely translates to solving Chapman-Kolmogorov
transit integral in downward direction.

Notes
- State machine function nr. 8a
- Includes initialization routines.
- TODO: Make multi-core
"""
function attemptCliqInitDown_StateMachine(csmc::CliqStateMachineContainer)
  #
  # should never happen to
  setCliqDrawColor(csmc.cliq, "green")
  csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing

  # initialize clique in downward direction
  # not if parent also needs downward init message
  infocsm(csmc, "8a, needs down message -- attempt down init")
  prnt = getParent(csmc.tree, csmc.cliq)[1]
  dwinmsgs = prepCliqInitMsgsDown!(csmc.cliqSubFg, csmc.tree, prnt)

  if length(dwinmsgs) == 0
    infocsm(csmc, "attemptCliqInitDown_StateMachine, no can do, must wait for siblings to update parent.")
    # go to 8c
    return waitChangeOnParentCondition_StateMachine
  end

  cliqst = doCliqInitDown!(csmc.cliqSubFg, csmc.cliq, dwinmsgs)
  # TODO: transfer values changed in the cliques should be transfered to the tree in proc 1 here.

  # TODO: maybe this should be here?
  setCliqStatus!(csmc.cliq, cliqst)

  # TODO move out
  children = getChildren(csmc.tree, csmc.cliq)
  if areCliqChildrenNeedDownMsg(children)
    # set messages if children :needdownmsg
    infocsm(csmc, "8a, doCliqInitDown! -- must set messages for future down init")
    # construct init's up msg to place in parent from initialized separator variables
    msg = prepCliqInitMsgsUp(csmc.cliqSubFg, csmc.cliq) # , tree,

    infocsm(csmc, "8a, putting fake upinitmsg in this cliq, msgs labels $(collect(keys(msg)))")
    # set fake up and notify down status
    setCliqUpInitMsgs!(csmc.cliq, csmc.cliq.index, msg)
    # setCliqStatus!(csmc.cliq, cliqst)
    setCliqDrawColor(csmc.cliq, "sienna")
    csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing

    notifyCliqDownInitStatus!(csmc.cliq, cliqst)

    infocsm(csmc, "8a, near-end down init attempt, $cliqst.")
  end

  return attemptCliqInitUp_StateMachine
end


"""
    $SIGNATURES

Delay loop if waiting on upsolves to complete.

Notes
- State machine 7b
"""
function slowCliqIfChildrenNotUpsolved_StateMachine(csmc::CliqStateMachineContainer)

  # special case short cut
  cliqst = getCliqStatus(csmc.cliq)
  if cliqst == :needdownmsg
    infocsm(csmc, "7b, shortcut on cliq is a needdownmsg status.")
    return isCliqNull_StateMachine
  end
  childs = getChildren(csmc.tree, csmc.cliq)
  len = length(childs)
  @inbounds for i in 1:len
    if getCliqStatus(childs[i]) != :upsolved
      infocsm(csmc, "7b, wait condition on upsolve, i=$i.")
      wait(getSolveCondition(childs[i]))
      break
    end
  end

  # go to 4
  return isCliqNull_StateMachine
end

"""
    $SIGNATURES

Notes
- State machine function nr. 7
"""
function determineCliqNeedDownMsg_StateMachine(csmc::CliqStateMachineContainer)

  # fetch children status
  stdict = blockCliqUntilChildrenHaveUpStatus(csmc.tree, csmc.cliq)

  # hard assumption here on upsolve from leaves to root
  proceed = true
  # fetch status from children (should already be available -- i.e. should not block)
  for (clid, clst) in stdict
    infocsm(csmc, "7, check stdict children: clid=$(clid), clst=$(clst)")
    # :needdownmsg # 'send' downward init msg direction
    !(clst in [:initialized;:upsolved;:marginalized;:downsolved]) ? (proceed = false) : nothing
  end
  infocsm(csmc, "7, proceed=$(proceed), forceproceed=$(csmc.forceproceed)")


  if proceed || csmc.forceproceed
    csmc.forceproceed = false
    # return doCliqInferAttempt_StateMachine
    cliqst = getCliqStatus(csmc.cliq)
    infocsm(csmc, "7, status=$(cliqst), before attemptCliqInitDown_StateMachine")
    # d1,d2,cliqst = doCliqInitUpOrDown!(csmc.cliqSubFg, csmc.tree, csmc.cliq, isprntnddw)
    if cliqst == :needdownmsg && !isCliqParentNeedDownMsg(csmc.tree, csmc.cliq)
      return attemptCliqInitDown_StateMachine
    end
    return attemptCliqInitUp_StateMachine
  else
    return slowCliqIfChildrenNotUpsolved_StateMachine
  end
end

"""
    $SIGNATURES

Notes
- State machine function nr. 6c
"""
function blockCliqSiblingsParentChildrenNeedDown_StateMachine(csmc::CliqStateMachineContainer)
  # add blocking case when all siblings and parent :needdownmsg -- until parent :initialized
  infocsm(csmc, "7, check/block sibl&prnt :needdownmsg")
  blockCliqSiblingsParentNeedDown(csmc.tree, csmc.cliq)

  return determineCliqNeedDownMsg_StateMachine
end


"""
    $SIGNATURES

Notes
- State machine function nr. 5
"""
function blockUntilSiblingsStatus_StateMachine(csmc::CliqStateMachineContainer)
  infocsm(csmc, "5, blocking on parent until all sibling cliques have valid status")
  setCliqDrawColor(csmc.cliq, "turquoise")
  csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing

  cliqst = getCliqStatus(csmc.cliq)
  @info "$(current_task()) Clique $(csmc.cliq.index), 5, block on siblings cliq status=$(cliqst)"
  prnt = getParent(csmc.tree, csmc.cliq)
  infocsm(csmc, "5, block on siblings cliq")
  if length(prnt) > 0
    infocsm(csmc, "5, has parent clique=$(prnt[1].index)")
    blockCliqUntilChildrenHaveUpStatus(csmc.tree, prnt[1])
  end

  # go to 6c
  return blockCliqSiblingsParentChildrenNeedDown_StateMachine
end



"""
    $SIGNATURES

Notes
- State machine function nr. 4
"""
function isCliqNull_StateMachine(csmc::CliqStateMachineContainer)

  #must happen before if :null
  stdict = blockCliqUntilChildrenHaveUpStatus(csmc.tree, csmc.cliq)
  csmc.forceproceed = false

  prnt = getParent(csmc.tree, csmc.cliq)
  if 0 == length(prnt)
    # go to 7
    return determineCliqNeedDownMsg_StateMachine
  end

  return doesCliqNeeddownmsg_StateMachine
end


"""
    $SIGNATURES

Determine if any down messages are required.

Notes
- State machine function nr.4b
"""
function doesCliqNeeddownmsg_StateMachine(csmc::CliqStateMachineContainer)
  # csmc.forceproceed = false
  cliqst = getCliqStatus(csmc.cliq)

  # infocsm(csmc, "4b, get parent")
  # # get parent cliq
  # prnt = getParent(csmc.tree, csmc.cliq)

  if cliqst != :null
    if cliqst != :needdownmsg
      return blockCliqSiblingsParentChildrenNeedDown_StateMachine
    end
    # ELSE: skip by to what used to be 4b
  else
    # fetch (should not block)
    stdict = blockCliqUntilChildrenHaveUpStatus(csmc.tree, csmc.cliq)
    chstatus = collect(values(stdict))
    len = length(chstatus)

    # if all children needdownmsg
    if len > 0 && sum(chstatus .== :needdownmsg) == len
      # TODO maybe can happen where some children need more information?
      infocsm(csmc, "4, escalating to :needdownmsg since all children :needdownmsg")
      notifyCliqUpInitStatus!(csmc.cliq, :needdownmsg)
      setCliqDrawColor(csmc.cliq, "green")
      csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing

      # go to 5
      return blockUntilSiblingsStatus_StateMachine
    end

    # go to 6c
    return blockCliqSiblingsParentChildrenNeedDown_StateMachine
  end # != :null

  # if cliqst == :needdownmsg
    if areCliqChildrenNeedDownMsg(csmc.tree, csmc.cliq)
      infocsm(csmc, "4, must deal with child :needdownmsg")
      csmc.forceproceed = true
    else
      # go to 5
      return blockUntilSiblingsStatus_StateMachine
    end
  # end

  # go to 6c
  return blockCliqSiblingsParentChildrenNeedDown_StateMachine
end


"""
    $SIGNATURES

Build a sub factor graph for clique variables from the larger factor graph.

Notes
- State machine function nr.2
"""
function buildCliqSubgraph_StateMachine(csmc::CliqStateMachineContainer)
  # build a local subgraph for inference operations
  infocsm(csmc, "2, build subgraph")
  syms = getCliqAllVarSyms(csmc.dfg, csmc.cliq)
  csmc.cliqSubFg = buildSubgraphFromLabels(csmc.dfg, syms)
  return isCliqNull_StateMachine
end


"""
    $SIGNATURES

Notify possible parent if clique is upsolved and exit the state machine.

Notes
- State machine function nr.1
"""
function isCliqUpSolved_StateMachine(csmc::CliqStateMachineContainer)
  infocsm(csmc, "1, isCliqUpSolved_StateMachine")
  cliqst = getCliqStatus(csmc.cliq)
  # lbl = cliq.attributes["label"]

  if csmc.incremental && cliqst in [:upsolved; :downsolved; :marginalized]
    # prep and send upward message
    prnt = getParent(csmc.tree, csmc.cliq)
    if length(prnt) > 0
      # not a root clique
      # construct init's up msg to place in parent from initialized separator variables
      msg = prepCliqInitMsgsUp(csmc.dfg, csmc.tree, csmc.cliq)
      setCliqUpInitMsgs!(prnt[1], csmc.cliq.index, msg)
      notifyCliqUpInitStatus!(csmc.cliq, cliqst)
    end
    return IncrementalInference.exitStateMachine
  end
  return buildCliqSubgraph_StateMachine
end


"""
    $SIGNATURES

EXPERIMENTAL: perform upward inference using a state machine solution approach.

Notes:
- will call on values from children or parent cliques
- can be called multiple times
- Assumes all cliques in tree are being solved simultaneously and in similar manner.
- State machine rev.1 -- copied from first TreeBasedInitialization.jl.
- Doesn't do partial initialized state properly yet.
"""
function cliqInitSolveUpByStateMachine!(dfg::G,
                                        tree::BayesTree,
                                        cliq::Graphs.ExVertex;
                                        drawtree::Bool=false,
                                        show::Bool=false,
                                        incremental::Bool=true,
                                        limititers::Int=-1,
                                        recordhistory::Bool=false  ) where G <: AbstractDFG
  #
  children = Graphs.ExVertex[]
  for ch in Graphs.out_neighbors(cliq, tree.bt)
    push!(children, ch)
  end
  prnt = getParent(tree, cliq)

  csmc = CliqStateMachineContainer(dfg, initfg(), tree, cliq, prnt, children, false, incremental, drawtree)

  statemachine = StateMachine{CliqStateMachineContainer}(next=isCliqUpSolved_StateMachine)
  while statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory); end
  statemachine.history
end



#
