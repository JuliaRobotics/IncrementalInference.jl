

"""
    $SIGNATURES

Return dict of all histories in a Bayes Tree.
"""
function getTreeCliqsSolverHistories(fg::G,
                                     tree::BayesTree)::Dict{Symbol, CSMHistory} where G <: AbstractDFG
  #
  fsy = getTreeAllFrontalSyms(fg, tree)
  histories = Dict{Symbol, CSMHistory}()
  for fs in fsy
    hist = getCliqSolveHistory(tree, fs)
    if length(hist) > 0
      histories[fs] = hist
    end
  end
  return histories
end


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
function printCliqHistorySummary(fid, hist::Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}})
  if length(hist) == 0
    @warn "printCliqHistorySummary -- No CSM history found."
  end
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
    println(fid, first)
  end
  nothing
end

function printCliqHistorySummary(hist::Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}})
  printCliqHistorySummary(stdout, hist)
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
                                step::Int)
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
function animateCliqStateMachines(tree::BayesTree,
                                  cliqsyms::Vector{Symbol};
                                  frames::Int=100  )
  #
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
    retval = animateStateMachineHistoryByTime(hist, frames=frames, folder="caesar/animatecsm/cliq$sym", title="$sym", startT=startT, stopT=stopT, rmfirst=false)
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
function solveCliqWithStateMachine!(dfg::G,
                                    tree::BayesTree,
                                    frontal::Symbol;
                                    iters::Int=200,
                                    downsolve::Bool=true,
                                    recordhistory::Bool=false,
                                    verbose::Bool=false,
                                    nextfnc::Function=isCliqUpSolved_StateMachine,
                                    prevcsmc::Union{Nothing,CliqStateMachineContainer}=nothing) where G <: AbstractDFG
  #
  cliq = whichCliq(tree, frontal)
  children = Graphs.ExVertex[]
  for ch in Graphs.out_neighbors(cliq, tree.bt)
    push!(children, ch)
  end
  prnt = getParent(tree, cliq)
  csmc = isa(prevcsmc, Nothing) ? CliqStateMachineContainer(dfg, initfg(), tree, cliq, prnt, children, false, true, true, downsolve) : prevcsmc
  statemachine = StateMachine{CliqStateMachineContainer}(next=nextfnc)
  while statemachine(csmc, verbose=verbose, iterlimit=iters, recordhistory=recordhistory); end
  statemachine, csmc
end

"""
    $SIGNATURES

Animate multiple clique state machines on the same graphviz visualization.  Renders according to
linear time for all provided histories.

Example:
```julia
using Caesar

# build a factor graph
fg = initfg()
# addVariable!(...)
# addFactor!(...)
# ...

fsy = getTreeAllFrontalSyms(fg, tree) # for later use
# perform inference to find the factor graph marginal posterior estimates
tree, smt, hist = solveTree!(fg, recordcliqs=fsy)

# generate frames in standard location /tmp/caesar/csmCompound/
#  requires: sudo apt-get install graphviz
csmAnimate(fg, tree, fsy, frames=500)

# to render and show from default location (requires )
#  requires: sudo apt-get install ffmpeg vlc
Base.rm("/tmp/caesar/csmCompound/out.mp4")
run(`ffmpeg -r 10 -i /tmp/caesar/csmCompound/csm_%d.png -c:v libx264 -vf fps=25 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" /tmp/caesar/csmCompound/out.mp4`)
run(`vlc /tmp/caesar/csmCompound/out.mp4`)
```
"""
function csmAnimate(fg::G,
                    tree::BayesTree,
                    cliqsyms::Vector{Symbol};
                    frames::Int=100,
                    rmfirst::Bool=true  ) where G <: AbstractDFG
  #
  hists = getTreeCliqsSolverHistories(fg,tree)

  startT = Dates.now()
  stopT = Dates.now()

  # get start and stop times across all cliques
  first = true
  for (csym, hist) in hists
    # global startT, stopT
    @show csym
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
  if rmfirst
    @warn "Removing /tmp/caesar/csmCompound/ in preparation for new frames."
    Base.rm("/tmp/caesar/csmCompound/", recursive=true, force=true)
  end
  animateStateMachineHistoryByTimeCompound(hists, startT, stopT, folder="caesar/csmCompound", frames=500)
end

"""
    $SIGNATURES

Return true if both, i.) this clique requires more downward information, ii.) more
downward message information could potentially become available.

Notes
- Delay initialization to the last possible moment.

Dev Notes:

Determine clique truely isn't able to proceed any further:
- should be as self reliant as possible (using clique's status as indicator)
- change status to :mustinitdown if have only partial beliefs so far:
  - combination of status, while partials belief siblings are not :mustinitdown
"""
function getCliqSiblingsPartialNeeds(tree::BayesTree, cliq::Graphs.ExVertex, prnt, dwinmsgs::Dict)
  # which incoming messages are partials
  hasPartials = Dict{Symbol, Int}()
  for (sym, tmsg) in dwinmsgs
    # assuming any down message per label that is not partial voids further partial consideration
    if sum(tmsg[2]) > 0
      if !haskey(hasPartials, sym)
        hasPartials[sym] = 0
      end
      hasPartials[sym] += 1
    end
  end
  partialKeys = collect(keys(hasPartials))

  ## determine who might be able to help init this cliq
  # check sibling separator sets against this clique's separator
  sibs = getCliqSiblings(tree, cliq)

  @info "getCliqSiblingsPartialNeeds -- CHECK PARTIAL"
  # identify which cliques might have useful information
  localsep = getCliqSeparatorVarIds(cliq)
  seps = Dict{Int, Vector{Symbol}}()
  for si in sibs
    # @show si.attributes["label"]
    mighthave = intersect(getCliqSeparatorVarIds(si), localsep)
    if length(mighthave) > 0
      seps[si.index] = mighthave
      # @show getCliqStatus(si)
      if getCliqStatus(si) in [:initialized; :null; :needdownmsg]
        # partials treated special -- this is slightly hacky
        # @show partialKeys, mighthave
        if length(intersect(localsep, partialKeys)) > 0 && length(mighthave) > 0
          # this sibling might have info to delay about
          setCliqDrawColor(cliq,"magenta")
          return true
        end
      end
    end
  end
  # determine if those cliques will / or will not be able to provide more info
  # when does clique change to :mustinitdown

  # TODO Should priorize solve order for mustinitdown with lowest dependency first

  # default
  return false
end


  # determine if any siblings might still hold promise
  # candidates = 0
  # for idx in collect(keys(seps))
  #   sibcliq = tree.cliques[idx]
  #   if getCliqStatus(sibcliq) in [:initialized; :null; :needdownmsg]
  #     @show sibcliq.attributes["label"],
  #     candidates += 1
  #   end
  # end
