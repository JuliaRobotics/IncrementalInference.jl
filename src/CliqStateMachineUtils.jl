
"""
    $SIGNATURES

Return clique state machine history from `tree` if it was solved with `recordcliqs`.

Notes
- Cliques are identified by front variable `::Symbol` which are always unique across the cliques.
"""
function getCliqSolveHistory(cliq::Graphs.ExVertex)
  getData(cliq).statehistory
end
function getCliqSolveHistory(tree::BayesTree, frntal::Symbol)
  cliq = whichCliq(tree, frntal)
  getCliqSolveHistory(cliq)
end

"""
    $SIGNATURES

Print a short summary of state machine history for a clique solve.

Related:

getTreeAllFrontalSyms, getCliqSolveHistory, animateCliqStateMachines
"""
function printCliqHistorySummary(hist::Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}})
  for hi in hist
    first = (split(string(hi[1]), 'T')[end])*" "
    len = length(first)
    for i in len:13  first = first*" "; end
    first = first*string(hi[2])
    len = length(first)
    for i in len:17  first = first*" "; end
    first = first*string(getCliqStatus(hi[4].cliq))
    len = length(first)
    for i in len:30  first = first*" "; end
    nextfn = split(split(string(hi[3]),'.')[end], '_')[1]
    lenf = length(nextfn)
    nextfn = 20 < lenf ? nextfn[1:20]*"." : nextfn
    first = first*nextfn
    len = length(first)
    for i in len:52  first = first*" "; end
    first = first*string(hi[4].forceproceed)
    len = length(first)
    for i in len:58  first = first*" "; end
    if 0 < length(hi[4].parentCliq)
      first = first*string(getCliqStatus(hi[4].parentCliq[1]))
    else
      first = first*"----"
    end
    first = first*" | "
    if 0 < length(hi[4].childCliqs)
      for ch in hi[4].childCliqs
        first = first*string(getCliqStatus(ch))*" "
      end
    end
    println(first)
  end
  nothing
end

function printCliqHistorySummary(cliq::Graphs.ExVertex)
  hist = getCliqSolveHistory(cliq)
  printCliqHistorySummary(hist)
end

function printCliqHistorySummary(tree::BayesTree, frontal::Symbol)
  hist = getCliqSolveHistory(tree, frontal)
  printCliqHistorySummary(hist)
end


"""
  $SIGNATURES

Repeat a solver state machine step without changing history or primary values.

printCliqHistorySummary, getCliqSolveHistory, cliqHistFilterTransitions
"""
function sandboxCliqResolveStep(tree::BayesTree,
                                frontal::Symbol,
                                step::Int  )
  #
  hist = getCliqSolveHistory(tree, frontal)
  return sandboxStateMachineStep(hist, step)
end




"""
    $SIGNATURES

Draw many images in '/tmp/?/csm_%d.png' representing time synchronized state machine
events for cliques `cliqsyms::Vector{Symbol}`.

Notes
- State history must have previously been recorded (stored in tree cliques).

Related

printCliqHistorySummary
"""
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

"""
    $SIGNATURES

Return state machine transition steps from history such that the `nextfnc::Function`.

Related:

getCliqSolveHistory, printCliqHistorySummary, filterHistAllToArray, sandboxCliqResolveStep
"""
function cliqHistFilterTransitions(hist::Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}, nextfnc::Function)
  ret = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}()
  for hi in hist
    if hi[3] == nextfnc
      push!(ret, hi)
    end
  end
  return ret
end

"""
    $SIGNATURES

Return state machine transition steps from all cliq histories with transition `nextfnc::Function`.

Related:

getCliqSolveHistory, printCliqHistorySummary, cliqHistFilterTransitions, sandboxCliqResolveStep
"""
function filterHistAllToArray(tree::BayesTree, frontals::Vector{Symbol}, nextfnc::Function)
  ret = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}()
  for sym in frontals
    hist = getCliqSolveHistory(tree, sym)
    fih = cliqHistFilterTransitions(hist, nextfnc)
    for fi in fih
      push!(ret, fi)
    end
   end
  return ret
end

"""
    $SIGNATURES

Standalone state machine solution for a single clique.

Related:

initInferTreeUp!
"""
function solveCliqWithStateMachine!(fg::FactorGraph,
                                    tree::BayesTree,
                                    frontal::Symbol;
                                    iters::Int=200,
                                    recordhistory::Bool=false,
                                    verbose::Bool=false  )
  #
  cliq = whichCliq(tree, frontal)
  children = Graphs.ExVertex[]
  for ch in Graphs.out_neighbors(cliq, tree.bt)
    push!(children, ch)
  end
  prnt = getParent(tree, cliq)
  csmc = CliqStateMachineContainer(fg, initfg(), tree, cliq, prnt, children, false, true, true)
  statemachine = StateMachine{CliqStateMachineContainer}(next=isCliqUpSolved_StateMachine)
  while statemachine(csmc, verbose=verbose, iterlimit=iters, recordhistory=recordhistory); end
  statemachine, csmc
end
