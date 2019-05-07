# clique state machine for tree based initialization and inference

mutable struct CliqStateMachineContainer
  fg::FactorGraph
  tree::BayesTree
  cliq::Graphs.ExVertex
  cliqSubFg::FactorGraph
  # TODO: bad flags that must be removed
  proceed::Bool
  forceproceed::Bool
  tryonce::Bool
  incremental::Bool
  drawtree::Bool
end


"""
    $SIGNATURES

Notes
- State machine function nr.9
"""
function finishCliqSolveCheck_StateMachine(csmc::CliqStateMachineContainer)
  cliqst = getCliqStatus(csmc.cliq)
  if cliqst == :upsolved
    @info "$(current_task()) Clique $(csmc.cliq.index), going for transferUpdateSubGraph!"
    frsyms = Symbol[getSym(csmc.cliqSubFg, varid) for varid in getCliqFrontalVarIds(csmc.cliq)]
    transferUpdateSubGraph!(csmc.fg, csmc.cliqSubFg, frsyms)
    return IncrementalInference.exitStateMachine
  elseif cliqst == :initialized
    # @info "$(current_task()) Clique $(cliq.index), set update down init messages: "  # OBSOLETE
    setCliqDrawColor(csmc.cliq, "sienna")
  else
    @info "$(current_task()) Clique $(csmc.cliq.index), init not complete and should wait on init down message."
    setCliqDrawColor(csmc.cliq, "green")
    csmc.tryonce = true # TODO, potential problem with trying to downsolve
  end
  csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing

  return whileCliqNotSolved_StateMachine
end


"""
    $SIGNATURES

Notes
- State machine function nr. 8
"""
function doCliqInferAttempt_StateMachine(csmc::CliqStateMachineContainer)
  setCliqDrawColor(csmc.cliq, "red")
  cliqst = getCliqStatus(csmc.cliq)
  csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing
  # evaluate according to cliq status
  isprntnddw = isCliqParentNeedDownMsg(csmc.tree, csmc.cliq)
  @info "$(current_task()) Clique $(csmc.cliq.index), proceed: $(cliqst), isCliqParentNeedDownMsg(tree, cliq)=$(isprntnddw), areCliqChildrenNeedDownMsg(tree, cliq)=$(areCliqChildrenNeedDownMsg(csmc.tree, csmc.cliq))"
  d1,d2,cliqst = doCliqInitUpOrDown!(csmc.cliqSubFg, csmc.tree, csmc.cliq, isprntnddw)

  return finishCliqSolveCheck_StateMachine
end


"""
    $SIGNATURES

Notes
- State machine function nr. 7
"""
function determineCliqNeedDownMsg_StateMachine(csmc::CliqStateMachineContainer)

  cliqst = getCliqStatus(csmc.cliq)
  lbl = csmc.cliq.attributes["label"]
  stdict = Dict{Int, Symbol}()

  # promote if longer down chain of :needdownmsg
  if cliqst == :null
    @info "$(current_task()) Clique $(csmc.cliq.index), determineCliqNeedDownMsg -- blocking until child cliques have status, cliqst=$(cliqst)"
    stdict = blockCliqUntilChildrenHaveUpStatus(csmc.tree, csmc.cliq)
    # TODO stdict here is just to get the status of child cliques
    @info "$(current_task()) Clique $(csmc.cliq.index) continue, children all have status"

    chstatus = collect(values(stdict))
    len = length(chstatus)
    if len > 0 && sum(chstatus .== :needdownmsg) == len
      # TODO maybe can happen where some children need more information?
      @info "$(current_task()) Clique $(csmc.cliq.index) | $lbl | escalating to :needdownmsg since all children :needdownmsg"
      notifyCliqUpInitStatus!(csmc.cliq, :needdownmsg)
      # setCliqStatus!(cliq, :needdownmsg)
      cliqst = getCliqStatus(csmc.cliq) ## TODO: likely not required since cliqst already exists
      setCliqDrawColor(csmc.cliq, "green")
      csmc.tryonce = true
    end

    # wait if child branches still solving -- must eventually upsolve this clique
    if len > 0 && sum(chstatus .!= :upsolved) > 0
      @info "$(current_task()) Clique $(csmc.cliq.index) | $lbl | sleeping until all children finish upward inference"
      sleep(0.1)
    end
  end
  # hard assumption here on upsolve from leaves to root
  csmc.proceed = true

  # TODO not sure if we want stdict from cliq or prnt???
  for (clid, clst) in stdict
    @info "$(current_task()) Clique $(csmc.cliq.index), check stdict: clid=$(clid), clst=$(clst)"
    # :needdownmsg # 'send' downward init msg direction
    # :initialized # @warn "something might not be right with init of clid=$clid"
    !(clst in [:initialized;:upsolved;:marginalized;:downsolved]) ? (csmc.proceed = false) : nothing
  end
  @info "$(current_task()) Clique $(csmc.cliq.index), proceed=$(csmc.proceed), tryonce=$(csmc.tryonce), clst=$(cliqst)"

  # add blocking case when all siblings and parent :needdownmsg -- until parent :initialized
  @info "$(current_task()) Clique $(csmc.cliq.index), check block if siblings & parent have :needdownmsg status? clst=$(cliqst), proceed=$(csmc.proceed), forceproceed=$(csmc.forceproceed)."
  blockCliqSiblingsParentNeedDown(csmc.tree, csmc.cliq)

  # add case for if children are blocked on need down msg
  if getCliqStatus(csmc.cliq) == :initialized && areCliqChildrenNeedDownMsg(csmc.tree, csmc.cliq)
    sleep(0.1)
  end

  if csmc.proceed || csmc.forceproceed
    return doCliqInferAttempt_StateMachine
  else
    return whileCliqNotSolved_StateMachine
  end
end


"""
    $SIGNATURES

Notes
- State machine function nr. 6
"""
function blockUntilChildrenStatus_StateMachine(csmc::CliqStateMachineContainer)
  cliqst = getCliqStatus(csmc.cliq)
  @info "$(current_task()) Clique $(csmc.cliq.index), blockUntilChildrenStatus_StateMachine -- blocking until child cliques have status, cliqst=$(cliqst)"
  blockCliqUntilChildrenHaveUpStatus(csmc.tree, csmc.cliq)
  @info "$(current_task()) Clique $(csmc.cliq.index) continue, children all have status"

  return determineCliqNeedDownMsg_StateMachine
end

"""
    $SIGNATURES

Notes
- State machine function nr. 5
"""
function blockUntilSiblingsStatus_StateMachine(csmc::CliqStateMachineContainer)
  prnt = getParent(csmc.tree, csmc.cliq)
  if length(prnt) > 0
    blockCliqUntilChildrenHaveUpStatus(csmc.tree, prnt[1])
  end
  return blockUntilChildrenStatus_StateMachine
end

"""
    $SIGNATURES

Determine if any down messages are required.

Notes
- State machine function nr.4
"""
function doesCliqNeeddownmsg_StateMachine(csmc::CliqStateMachineContainer)
  csmc.tryonce = false
  csmc.forceproceed = false
  cliqst = getCliqStatus(csmc.cliq)
  # stdictprnt = Dict{Int, Symbol}()
  # get parent cliq
  prnt = getParent(csmc.tree, csmc.cliq)

  # @info "$(current_task()) Clique $(cliq.index), status $cliqst -- top of while loop"
  if cliqst == :needdownmsg && length(prnt) > 0
    # wait here until all children have a valid status
    if !areCliqChildrenNeedDownMsg(csmc.tree, csmc.cliq)
      @info "$(current_task()) Clique $(csmc.cliq.index), blocking on parent until all sibling cliques have valid status"
      setCliqDrawColor(csmc.cliq, "turquoise")
      csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing
      return blockUntilSiblingsStatus_StateMachine
      # stdictprnt = blockCliqUntilChildrenHaveUpStatus(tree, prnt[1])
    else
      @warn "$(current_task()) Clique $(csmc.cliq.index), WIP must deal with child :needdownmsg"
      csmc.forceproceed = true
    end
  end
  return blockUntilChildrenStatus_StateMachine
end

"""
    $SIGNATURES

Determine if necessary to continue with solve attempt of this `csmc.cliq`.

Notes
- State machine function nr.3
- TODO: LIKELY MISSING A NOTIFY STEP -- should probably point at `isCliqUpSolved_StateMachine`.
"""
function whileCliqNotSolved_StateMachine(csmc::CliqStateMachineContainer)
  cliqst = getCliqStatus(csmc.cliq)
  if csmc.tryonce || !(cliqst in [:upsolved; :downsolved; :marginalized])
    return doesCliqNeeddownmsg_StateMachine
  else
    @info "Exit cliq state machine at whileCliqNotSolved_StateMachine"
    return IncrementalInference.exitStateMachine
  end
end

"""
    $SIGNATURES

Build a sub factor graph for clique variables from the larger factor graph.

Notes
- State machine function nr.2
"""
function buildCliqSubgraph_StateMachine(csmc::CliqStateMachineContainer)
  # build a local subgraph for inference operations
  syms = getCliqAllVarSyms(csmc.fg, csmc.cliq)
  csmc.cliqSubFg = buildSubgraphFromLabels(csmc.fg, syms)
  return whileCliqNotSolved_StateMachine
end


"""
    $SIGNATURES

Notify possible parent if clique is upsolved and exit the state machine.

Notes
- State machine function nr.1
"""
function isCliqUpSolved_StateMachine(csmc::CliqStateMachineContainer)
  @info "$(current_task()) Clique $(csmc.cliq.index), isCliqUpSolved_StateMachine"
  cliqst = getCliqStatus(csmc.cliq)
  # lbl = cliq.attributes["label"]

  if csmc.incremental && cliqst in [:upsolved; :downsolved; :marginalized]
    # prep and send upward message
    prnt = getParent(csmc.tree, csmc.cliq)
    if length(prnt) > 0
      # not a root clique
      # construct init's up msg to place in parent from initialized separator variables
      msg = prepCliqInitMsgsUp!(csmc.fg, csmc.tree, csmc.cliq)
      setCliqUpInitMsgs!(prnt[1], csmc.cliq.index, msg)
      notifyCliqUpInitStatus!(csmc.cliq, cliqst)
      # @info "$(current_task()) Clique $(csmc.cliq.index), skip computation on status=$(cliqst), but did prepare/notify upward message"
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
function cliqInitSolveUpByStateMachine!(fg::FactorGraph,
                                        tree::BayesTree,
                                        cliq::Graphs.ExVertex;
                                        drawtree::Bool=false,
                                        show::Bool=false,
                                        incremental::Bool=true,
                                        limititers::Int=-1,
                                        recordhistory::Bool=false  )
  #
  csmc = CliqStateMachineContainer(fg, tree, cliq, initfg(), true, false, false, incremental, drawtree)

  statemachine = StateMachine{CliqStateMachineContainer}(next=isCliqUpSolved_StateMachine)
  while statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory); end
  statemachine.history
end






#
