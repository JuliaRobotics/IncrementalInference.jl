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

  with_logger(csmc.logger) do
    @info "$tmt | $(current_task()) cliq $(csmc.cliq.index), $lbl1, $(cliqst) -- "*str
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
  csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing

  # get down msg from parent
  prnt = getParent(csmc.tree, csmc.cliq)
  dwnmsgs = getDwnMsgs(prnt[1])
  multiproc = true

  # call down inference, TODO multiproc
  if multiproc
    @info "GOING MULTIPROC DWN"
    cliqc = deepcopy(csmc.cliq)
    cliqcd = getData(cliqc)
    # redirect to new unused so that CAN be serialized
    cliqcd.initUpChannel = Channel{Symbol}(1)
    cliqcd.initDownChannel = Channel{Symbol}(1)
    cliqcd.solveCondition = Condition()
    cliqcd.statehistory = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}()
    drt = remotecall_fetch(downGibbsCliqueDensity, upp2(), csmc.cliqSubFg, cliqc, dwnmsgs, 100, 3, false)
  else
    drt = downGibbsCliqueDensity(csmc.cliqSubFg, csmc.cliq, dwnmsgs, 100, 3, false)
  end
  csmc.dodownsolve = false

  # update clique with new status
  updateFGBT!(csmc.cliqSubFg, csmc.tree, csmc.cliq.index, drt, dbg=false, fillcolor="lightblue")
  setCliqStatus!(csmc.cliq, :downsolved) # should be a notify
  notifyCliqDownInitStatus!(csmc.cliq, :downsolved)

  csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing
  # and finished
  return IncrementalInference.exitStateMachine
end


"""
    $SIGNATURES

Direct state machine to continue with downward solve or exit.

Notes
- State machine function nr. 10
"""
function determineCliqIfDownSolve_StateMachine(csmc::CliqStateMachineContainer)
    infocsm(csmc, "10, determineCliqIfDownSolve_StateMachine, csmc.dodownsolve=$(csmc.dodownsolve).")
  # finished and exit downsolve
  if !csmc.dodownsolve
    infocsm(csmc, "10, determineCliqIfDownSolve_StateMachine -- shortcut exit since downsolve not required.")
    return IncrementalInference.exitStateMachine
  end

  # yes, continue with downsolve
  setCliqDrawColor(csmc.cliq, "turquoise")
  csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing

  # block here until parent is downsolved
  prnt = getParent(csmc.tree, csmc.cliq)
  if length(prnt) > 0
    blockCliqUntilParentDownSolved(prnt[1])
  else
    # this is the root clique, so assume already downsolved -- only special case
    setCliqDrawColor(csmc.cliq, "lightblue")
    dwnmsgs = getCliqDownMsgsAfterDownSolve(csmc.cliqSubFg, csmc.cliq)
    setDwnMsg!(csmc.cliq, dwnmsgs)
    setCliqStatus!(csmc.cliq, :downsolved)
	csmc.dodownsolve = false
    notifyCliqDownInitStatus!(csmc.cliq, :downsolved)
    return IncrementalInference.exitStateMachine
  end

  # go to 11
  return doCliqDownSolve_StateMachine
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
    transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frsyms) # TODO what about down solve??
    # go to 10
    return determineCliqIfDownSolve_StateMachine # IncrementalInference.exitStateMachine
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
    infocsm(csmc, "8c, waitChangeOnParentCondition_StateMachine, wait on parent for condition notify.")
    wait(getSolveCondition(prnt[1]))
  else
    infocsm(csmc, "8c, waitChangeOnParentCondition_StateMachine, cannot wait on parent for condition notify.")
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

  if csmc.delay
    infocsm(csmc, "8b, attemptCliqInitUp, delay required -- sleeping for 10s.")
    sleep(30)
  end

  cliqst = getCliqStatus(csmc.cliq)

  infocsm(csmc, "8b, attemptCliqInitUp, !areCliqChildrenNeedDownMsg()=$(!areCliqChildrenNeedDownMsg(csmc.tree, csmc.cliq))" )
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

  # determine if more info is needed for partial
  partialneedsmore = getCliqSiblingsPartialNeeds(csmc.tree, csmc.cliq, prnt, dwinmsgs)

  if length(dwinmsgs) == 0 || partialneedsmore
    infocsm(csmc, "8a, attemptCliqInitDown_StateMachine, no can do, must wait for siblings to update parent.")
    # go to 8c
    return waitChangeOnParentCondition_StateMachine
  end

  ## TODO deal with partial inits only, either delay or continue at end...
  # find intersect between downinitmsgs and local clique variables
  # if only partials available, then

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
    if !(getCliqStatus(childs[i]) in [:upsolved;:uprecycled;:marginalized])
      infocsm(csmc, "7b, wait condition on upsolve, i=$i, ch_lbl=$(getCliqFrontalVarIds(childs[i])[1]).")
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
    !(clst in [:initialized;:upsolved;:marginalized;:downsolved;:uprecycled]) ? (proceed = false) : nothing
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
  setCliqDrawColor(csmc.cliq, "darkviolet")
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

  prnt = getParent(csmc.tree, csmc.cliq)
  infocsm(csmc, "4, isCliqNull_StateMachine, csmc.incremental=$(csmc.incremental), len(prnt)=$(length(prnt))")
  #must happen before if :null
  stdict = blockCliqUntilChildrenHaveUpStatus(csmc.tree, csmc.cliq)
  csmc.forceproceed = false

  # for recycle computed clique values case
  if csmc.incremental && getCliqStatus(csmc.oldcliqdata) == :downsolved
    csmc.incremental = false
    # might be able to recycle the previous clique solve, go to 0b
    return checkChildrenAllUpRecycled_StateMachine
  end

  if 0 == length(prnt)
    # go to 7
    return determineCliqNeedDownMsg_StateMachine
  end

  # go to 4b
  return doesCliqNeeddownmsg_StateMachine
end


"""
    $SIGNATURES

Determine if any down messages are required.

Notes
- State machine function nr.4b
"""
function doesCliqNeeddownmsg_StateMachine(csmc::CliqStateMachineContainer)

  cliqst = getCliqStatus(csmc.cliq)
  infocsm(csmc, "4b, doesCliqNeeddownmsg_StateMachine, cliqst=$cliqst")

  if cliqst != :null
    if cliqst != :needdownmsg
      return blockCliqSiblingsParentChildrenNeedDown_StateMachine
    end
  else
    # fetch (should not block)
    stdict = blockCliqUntilChildrenHaveUpStatus(csmc.tree, csmc.cliq)
    chstatus = collect(values(stdict))
    len = length(chstatus)

    # if all children needdownmsg
    if len > 0 && sum(chstatus .== :needdownmsg) == len
      # TODO maybe can happen where some children need more information?
      infocsm(csmc, "4b, escalating to :needdownmsg since all children :needdownmsg")
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
  syms = getCliqAllVarSyms(csmc.dfg, csmc.cliq)
  infocsm(csmc, "2, build subgraph syms=$(syms)")
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

  if cliqst in [:upsolved; :downsolved; :marginalized; :uprecycled]  #moved to 4 --- csmc.incremental &&
    # prep and send upward message
    prnt = getParent(csmc.tree, csmc.cliq)
    if length(prnt) > 0
      # not a root clique
      # construct init's up msg to place in parent from initialized separator variables
      msg = prepCliqInitMsgsUp(csmc.dfg, csmc.tree, csmc.cliq)
      setCliqUpInitMsgs!(prnt[1], csmc.cliq.index, msg)
      notifyCliqUpInitStatus!(csmc.cliq, cliqst)
    end
    #go to 10
    return determineCliqIfDownSolve_StateMachine
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
    chst = getCliqStatus(ch)
    if chst in [:uprecycled; :marginalized]
      push!(count, 1)
    end
  end
  infocsm(csmc, "0b, checkChildrenAllUpRecycled_StateMachine -- length(chldr)=$(length(chldr)), sum(count)=$(sum(count))")

  # all children can be used for uprecycled -- i.e. no children have new information
  if sum(count) == length(chldr)
    # set up msg and exit go to 1
    setCliqStatus!(csmc.cliq, :uprecycled)
    setCliqDrawColor(csmc.cliq, "orange")
    csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing
    # go to 1
    return isCliqUpSolved_StateMachine
  end

  # return to regular solve, go to 2
  return buildCliqSubgraph_StateMachine
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
  # @show getCliqFrontalVarIds(csmc.oldcliqdata), getCliqStatus(csmc.oldcliqdata)
  infocsm(csmc, "0., $(csmc.incremental) ? :uprecycled => getCliqStatus(csmc.oldcliqdata)=$(getCliqStatus(csmc.oldcliqdata))")

  # check if should be trying and can recycle clique computations
  if csmc.incremental && getCliqStatus(csmc.oldcliqdata) == :downsolved
    # check if a subgraph will be needed later
    if csmc.dodownsolve
      # yes need subgraph and need more checks, so go to 2
      return buildCliqSubgraph_StateMachine
    else
       # one or two checks say yes, so go to 4
      return isCliqNull_StateMachine
    end
  end

  # nope, regular clique init-solve, go to 1
  return isCliqUpSolved_StateMachine
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
                                        N::Int=100,
										oldcliqdata::BayesTreeNodeData=emptyBTNodeData(),
                                        drawtree::Bool=false,
                                        show::Bool=false,
                                        incremental::Bool=true,
                                        limititers::Int=-1,
                                        downsolve::Bool=false,
                                        recordhistory::Bool=false,
                                        delay::Bool=false,
                                        logger::SimpleLogger=SimpleLogger(Base.stdout)) where {G <: AbstractDFG, AL <: AbstractLogger}
  #
  children = Graphs.ExVertex[]
  for ch in Graphs.out_neighbors(cliq, tree.bt)
    push!(children, ch)
  end
  prnt = getParent(tree, cliq)

  csmc = CliqStateMachineContainer(dfg, initfg(), tree, cliq, prnt, children, false, incremental, drawtree, downsolve, delay, oldcliqdata, logger)

  statemachine = StateMachine{CliqStateMachineContainer}(next=testCliqCanRecycled_StateMachine)
  while statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory); end
  statemachine.history
end



#
