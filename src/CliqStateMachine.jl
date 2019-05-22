# clique state machine for tree based initialization and inference


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
    frsyms = Symbol[getSym(csmc.cliqSubFg, varid) for varid in getCliqFrontalVarIds(csmc.cliq)]
    transferUpdateSubGraph!(csmc.fg, csmc.cliqSubFg, frsyms)
    return IncrementalInference.exitStateMachine
  elseif cliqst == :initialized
    setCliqDrawColor(csmc.cliq, "sienna")
    csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing
  else
    infocsm(csmc, "9, init not complete and should wait on init down message.")
    setCliqDrawColor(csmc.cliq, "green")
    # TODO, potential problem with trying to downsolve
    return doesCliqNeeddownmsg_StateMachine
  end

  return whileCliqNotSolved_StateMachine
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
  if cliqst in [:initialized; :null] && !areCliqChildrenNeedDownMsg(csmc.tree, csmc.cliq)
    setCliqDrawColor(csmc.cliq, "red")
    csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing
    cliqst = doCliqAutoInitUp!(csmc.cliqSubFg, csmc.tree, csmc.cliq)
  end

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
  setCliqDrawColor(csmc.cliq, "red")
  csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing

  # initialize clique in downward direction
  # not if parent also needs downward init message
  infocsm(csmc, "8a, needs down message -- attempt down init")
  prnt = getParent(csmc.tree, csmc.cliq)[1]
  dwinmsgs = prepCliqInitMsgsDown!(csmc.cliqSubFg, csmc.tree, prnt)

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

    infocsm(csmc, "8a, after down init attempt, $cliqst.")
  end

  # repeat the if a second time, is bad TODO
  if cliqst in [:initialized; :null] && !areCliqChildrenNeedDownMsg(csmc.tree, csmc.cliq)
    return attemptCliqInitUp_StateMachine
  end
  return finishCliqSolveCheck_StateMachine

  # return attemptCliqInitUp_StateMachine
end


"""
    $SIGNATURES

Do actual inference calculations, loosely translates to solving Chapman-Kolmogorov transit integral in
either up or downward direction, although some caveats on when which occurs.

Notes
- State machine function nr. 8
- Used both during downward initialization and upward initialization / full-solve.
"""
function doCliqInferAttempt_StateMachine(csmc::CliqStateMachineContainer)
  # visualization and debugging
  csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing

  # evaluate according to cliq status
  cliqst = getCliqStatus(csmc.cliq)

  infocsm(csmc, "8, status=$(cliqst), before attemptCliqInitDown_StateMachine")
  # d1,d2,cliqst = doCliqInitUpOrDown!(csmc.cliqSubFg, csmc.tree, csmc.cliq, isprntnddw)
  if cliqst == :needdownmsg && !isCliqParentNeedDownMsg(csmc.tree, csmc.cliq)
    return attemptCliqInitDown_StateMachine
  end

  # cliqst = getCliqStatus(csmc.cliq)
  if cliqst in [:initialized; :null] && !areCliqChildrenNeedDownMsg(csmc.tree, csmc.cliq)
    return attemptCliqInitUp_StateMachine
  end
  # either finished or going around again
  return finishCliqSolveCheck_StateMachine

  # return attemptCliqInitUp_StateMachine
end


function slowCliqIfChildrenNotUpsolved_StateMachine(csmc::CliqStateMachineContainer)
  childs = getChildren(csmc.tree, csmc.cliq)
  len = length(childs)
  tps = Vector{Bool}(len)
  fill!(tps, false)
  for i in 1:len
    ch = childs[i]
    if getCliqStatus(ch) == :upsolved
      tps[i] = true
    end
  end
  if sum(tps) != len
    infocsm(csmc, "7b, delay 0.2s since not all children are yet upsolved (hack).")
    sleep(0.2)
  end
  return whileCliqNotSolved_StateMachine
end

"""
    $SIGNATURES

Notes
- State machine function nr. 7
"""
function determineCliqNeedDownMsg_StateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "7, determineCliqNeedDownMsg_StateMachine start")
  cliqst = getCliqStatus(csmc.cliq)
  lbl = csmc.cliq.attributes["label"]
  stdict = Dict{Int, Symbol}()

  # promote if longer down chain of :needdownmsg
  if cliqst == :null
    infocsm(csmc, "7, determineCliqNeedDownMsg -- blocking until child cliques have status")
    stdict = blockCliqUntilChildrenHaveUpStatus(csmc.tree, csmc.cliq)
    # TODO stdict here is just to get the status of child cliques

    chstatus = collect(values(stdict))
    len = length(chstatus)
    if len > 0 && sum(chstatus .== :needdownmsg) == len
      # TODO maybe can happen where some children need more information?
      infocsm(csmc, "7, escalating to :needdownmsg since all children :needdownmsg")
      notifyCliqUpInitStatus!(csmc.cliq, :needdownmsg)
      # setCliqStatus!(cliq, :needdownmsg)
      cliqst = getCliqStatus(csmc.cliq) ## TODO: likely not required since cliqst already exists
      setCliqDrawColor(csmc.cliq, "green")
      csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing

      return blockUntilSiblingsStatus_StateMachine
    end

    # TODO: REMOVE
    # if len > 0 && sum(chstatus .!= :upsolved) > 0

  end

  # hard assumption here on upsolve from leaves to root
  proceed = true

  # TODO not sure if we want stdict from cliq or prnt???
  for (clid, clst) in stdict
    infocsm(csmc, "7, check stdict children: clid=$(clid), clst=$(clst)")
    # :needdownmsg # 'send' downward init msg direction
    !(clst in [:initialized;:upsolved;:marginalized;:downsolved]) ? (proceed = false) : nothing
  end
  infocsm(csmc, "7, proceed=$(proceed)")

  # add blocking case when all siblings and parent :needdownmsg -- until parent :initialized
  infocsm(csmc, "7, check/block sibl&prnt :needdownmsg, proceed=$(proceed), forceproceed=$(csmc.forceproceed).")
  blockCliqSiblingsParentNeedDown(csmc.tree, csmc.cliq)

  # # add case for if children are blocked on need down msg
  # # TODO: remove as soon possible
  # if getCliqStatus(csmc.cliq) == :initialized && areCliqChildrenNeedDownMsg(csmc.tree, csmc.cliq)
  #   sleep(0.1)
  # end

  if proceed || csmc.forceproceed
    return doCliqInferAttempt_StateMachine
  else
    return slowCliqIfChildrenNotUpsolved_StateMachine
  end
end


"""
    $SIGNATURES

Notes
- State machine function nr. 6
"""
function blockUntilChildrenStatus_StateMachine(csmc::CliqStateMachineContainer)
  cliqst = getCliqStatus(csmc.cliq)
  infocsm(csmc, "6, blockUntilChildrenStatus_StateMachine -- blocking until child cliques have status")
  blockCliqUntilChildrenHaveUpStatus(csmc.tree, csmc.cliq)
  infocsm(csmc, "6, continue, children all have status")

  return determineCliqNeedDownMsg_StateMachine
end

"""
    $SIGNATURES

Notes
- State machine function nr. 5
"""
function blockUntilSiblingsStatus_StateMachine(csmc::CliqStateMachineContainer)
  prnt = getParent(csmc.tree, csmc.cliq)
  infocsm(csmc, "5, block on siblings cliq")
  if length(prnt) > 0
    infocsm(csmc, "5, has parent clique=$(prnt[1].index)")
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

  infocsm(csmc, "4, get parent")
  # get parent cliq
  prnt = getParent(csmc.tree, csmc.cliq)

  # @info "$(current_task()) Clique $(cliq.index), status $cliqst -- top of while loop"
  if cliqst == :needdownmsg && length(prnt) > 0
    # wait here until all children have a valid status
    if !areCliqChildrenNeedDownMsg(csmc.tree, csmc.cliq)
      infocsm(csmc, "4, blocking on parent until all sibling cliques have valid status")
      setCliqDrawColor(csmc.cliq, "turquoise")
      csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing
      return blockUntilSiblingsStatus_StateMachine
    else
      infocsm(csmc, "4, WIP must deal with child :needdownmsg")
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
  infocsm(csmc, "3, whileCliqNotSolved_StateMachine")

  if !(cliqst in [:upsolved; :downsolved; :marginalized])  # csmc.tryonce ||
    return doesCliqNeeddownmsg_StateMachine
  else
    infocsm(csmc, "3, Exit cliq state machine at whileCliqNotSolved_StateMachine")
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
  infocsm(csmc, "2, build subgraph")
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
  infocsm(csmc, "1, isCliqUpSolved_StateMachine")
  cliqst = getCliqStatus(csmc.cliq)
  # lbl = cliq.attributes["label"]

  if csmc.incremental && cliqst in [:upsolved; :downsolved; :marginalized]
    # prep and send upward message
    prnt = getParent(csmc.tree, csmc.cliq)
    if length(prnt) > 0
      # not a root clique
      # construct init's up msg to place in parent from initialized separator variables
      msg = prepCliqInitMsgsUp(csmc.fg, csmc.tree, csmc.cliq)
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
  children = Graphs.ExVertex[]
  for ch in Graphs.out_neighbors(cliq, tree.bt)
    push!(children, ch)
  end
  prnt = getParent(tree, cliq)

  csmc = CliqStateMachineContainer(fg, initfg(), tree, cliq, prnt, children, false, false, incremental, drawtree)

  statemachine = StateMachine{CliqStateMachineContainer}(next=isCliqUpSolved_StateMachine)
  while statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory); end
  statemachine.history
end


"""
    $SIGNATURES

Return clique state machine history from `tree` if it was solved with `recordcliqs`.

Notes
- Cliques are identified by front variable `::Symbol` which are always unique across the cliques.
"""
function getCliqSolveHistory(tree::BayesTree, frntal::Symbol)
  cliq = whichCliq(tree, frntal)
  getData(cliq).statehistory
end

"""
    $SIGNATURES

Print a short summary of state machine history for a clique solve.
"""
function printCliqHistorySummary(tree::BayesTree, frontal::Symbol)
  hist = getCliqSolveHistory(tree, frontal)
  for hi in hist
    first = string(hi[1])
    len = length(first)
    for i in len:3  first = first*" "; end
    first = first*string(getCliqStatus(hi[3].cliq))
    len = length(first)
    for i in len:15  first = first*" "; end
    first = first*split(split(string(hi[2]),'.')[end], '_')[1]
    len = length(first)
    for i in len:40  first = first*" "; end
    first = first*string(hi[3].forceproceed)
    len = length(first)
    for i in len:46  first = first*" "; end
    first = first*string(hi[3].tryonce)
    println(first)
  end
  nothing
end


"""
  $SIGNATURES

Repeat a solver state machine step without changing history or primary values.
"""
function sandboxCliqResolveStep(tree::BayesTree,
                                frontal::Symbol,
                                step::Int  )
  #
  hist = getCliqSolveHistory(tree, frontal)
  return sandboxStateMachineStep(hist, step)
end





function animateCliqStateMachines(tree::BayesTree, cliqsyms::Vector{Symbol}; frames::Int=100)

  startT = Dates.now()
  stopT = Dates.now()

  # get start and stop times across all cliques
  first = true
  for sym in cliqsyms
    hist = getCliqSolveHistory(tree, sym)
    if hist[1][1] < startT
      startT = hist[1][1]
    end
    if first
      stopT = hist[end][1]
    end
    if stopT < hist[end][1]
      stopT= hist[end][1]
    end
  end

  # export all figures
  folders = String[]
  for sym in cliqsyms
    hist = getCliqSolveHistory(tree, sym)
    retval = animateStateMachineHistoryByTime(hist, frames=frames, folder="cliq$sym", title="$sym", startT=startT, stopT=stopT)
    push!(folders, "cliq$sym")
  end

  return folders
end



#
