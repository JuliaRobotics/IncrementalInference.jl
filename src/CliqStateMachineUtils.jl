

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

function printCliqHistorySummary(hists::Dict{Int,Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}, tree::BayesTree, sym::Symbol)
  hist = hists[getCliq(tree, sym).index]
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

printCliqSummary, printCliqHistorySummary, getCliqSolveHistory, cliqHistFilterTransitions
"""
function sandboxCliqResolveStep(tree::BayesTree,
                                frontal::Symbol,
                                step::Int)
  #
  hist = getCliqSolveHistory(tree, frontal)
    # clear Condition states to allow step solve
    cond = getSolveCondition(hist[step][4].cliq)
    cond.waitq = Any[]
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
      first = false
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
                                    nextfnc::Function=testCliqCanRecycled_StateMachine,
                                    prevcsmc::Union{Nothing,CliqStateMachineContainer}=nothing) where G <: AbstractDFG
  #
  cliq = whichCliq(tree, frontal)
  children = Graphs.ExVertex[]
  for ch in Graphs.out_neighbors(cliq, tree.bt)
    push!(children, ch)
  end
  prnt = getParent(tree, cliq)

  destType = (G <: InMemoryDFGTypes) ? G : InMemDFGType#GraphsDFG{SolverParams}

  csmc = isa(prevcsmc, Nothing) ? CliqStateMachineContainer(dfg, initfg(destType), tree, cliq, prnt, children, false, true, true, downsolve, false, getSolverParams(dfg)) : prevcsmc
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

# to render and show from default location (might require)
#  sudo apt-get install ffmpeg vlc

# .ogv [Totem Ubuntu default]
Base.rm("/tmp/caesar/csmCompound/out.ogv")
run(`ffmpeg -r 10 -i /tmp/caesar/csmCompound/csm_%d.png -c:v libtheora -vf fps=25 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -q 10 /tmp/caesar/csmCompound/out.ogv`)
run(`totem /tmp/caesar/csmCompound/out.ogv`)

# h.264 [VLC not default]
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
      first = false
    end
    if stopT < hist[end][1]
      stopT = hist[end][1]
    end
  end

  # export all figures
  if rmfirst
    @warn "Removing /tmp/caesar/csmCompound/ in preparation for new frames."
    Base.rm("/tmp/caesar/csmCompound/", recursive=true, force=true)
  end
  animateStateMachineHistoryByTimeCompound(hists, startT, stopT, folder="caesar/csmCompound", frames=frames)
end

"""
    $SIGNATURES

Convenience function to assign and make video of CSM state machine for `cliqs`.

Notes
- Probably several teething issues still (lower priority).
- Use `assignhist` if solver params async was true, or errored.

Related

csmAnimate, printCliqHistorySummary
"""
function makeCsmMovie(fg::G,
                      tree::BayesTree,
                      cliqs=ls(fg);
                      assignhist=nothing,
                      show::Bool=true,
                      filename::AS="/tmp/caesar/csmCompound/out.ogv",
                      frames::Int=1000 )::String where {G <: AbstractDFG, AS <: AbstractString}
  #
  if assignhist != nothing
    assignTreeHistory!(tree, assignhist)
  end
  csmAnimate(fg, tree, cliqs, frames=frames)
  # Base.rm("/tmp/caesar/csmCompound/out.ogv")
  run(`ffmpeg -r 10 -i /tmp/caesar/csmCompound/csm_%d.png -c:v libtheora -vf fps=25 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -q 10 $filename`)
  if show
    @async run(`totem $filename`)
  end
  filename
end




"""
    $SIGNATURES

Return true if this clique's down init should be delayed on account of prioritization among sibling separators.

Notes
- process described in issue #344

Dev Notes
- not priorizing order yet (TODO), just avoiding unsolvables at this time.
- Very closely related to getCliqSiblingsPartialNeeds -- refactor likely (NOTE).
- should precompute `allinters`.
"""
function getSiblingsDelayOrder(tree::BayesTree, cliq::Graphs.ExVertex, prnt, dwinmsgs::Dict; logger=ConsoleLogger())
  # when is a cliq upsolved
  solvedstats = Symbol[:upsolved; :marginalized; :uprecycled]

  # safety net double check
  cliqst = getCliqStatus(cliq)
  if cliqst in solvedstats
    with_logger(logger) do
      @warn "getSiblingsDelayOrder -- clique status should not be here with a solved cliqst=$cliqst"
    end
    return false
  end

  # get siblings separators
  sibs = getCliqSiblings(tree, cliq, true)
  ids = map(s->s.index, sibs)
  len = length(sibs)
  sibidx = collect(1:len)[ids .== cliq.index][1]
  seps = getCliqSeparatorVarIds.(sibs)
  lielbls = setdiff(ids, cliq.index)
  # get intersect matrix of siblings (should be exactly the same across siblings' csm)
  allinters = Array{Int,2}(undef, len, len)
  dwninters = Vector{Int}(undef, len)
  dwnkeys = collect(keys(dwinmsgs))
  with_logger(logger) do
    @info "getSiblingsDelayOrder -- number siblings=$(len), sibidx=$sibidx"
  end

  # sum matrix with all "up solved" rows and columns eliminated
  fill!(allinters, 0)
  for i in 1:len
    for j in i:len
      if i != j
        allinters[i,j] = length(intersect(seps[i],seps[j]))
      end
    end
    dwninters[i] = length(intersect(seps[i], dwnkeys))
  end

  # sum "across/over" rows, then columns (i.e. visa versa "along" columns then rows)
  rows = sum(allinters, dims=1)
  cols = sum(allinters, dims=2)

  with_logger(logger) do
      @info "getSiblingsDelayOrder -- allinters=$(allinters)"
      @info "getSiblingsDelayOrder -- rows=$(rows)"
      @info "getSiblingsDelayOrder -- rows=$(cols)"
  end

  # is this clique a non-zero row -- i.e. sum across columns? if not, no further special care needed
  if cols[sibidx] == 0
    with_logger(logger) do
      @info "getSiblingsDelayOrder -- cols[sibidx=$(sibidx))] == 0, no special care needed"
    end
    return false
  end

  # now determine if initializing from below or needdownmsg
  if cliqst in Symbol[:needdownmsg;]
    # be super careful about delay (true) vs pass (false) at this point -- might be partial too TODO
    # return true if delay beneficial to initialization accuracy

    # find which siblings this cliq epends on
    symm = allinters + allinters'
    maskcol = symm[:,sibidx] .> 0
    # lenm = length(maskcol)
    stat = Vector{Symbol}(undef, len)
    stillbusymask = fill(false, len)

    flush(logger.stream)

    # get each sibling status (entering atomic computation segment -- until wait command)
    stat .= getCliqStatus.(sibs) #[maskcol]

    ## (long down chain case)
    # need different behaviour when all remaining siblings are blocking with :needdownmsg
    remainingmask = stat .== :needdownmsg
    if sum(remainingmask) == length(stat)
      with_logger(logger) do
        @info "getSiblingsDelayOrder -- all blocking: sum(remainingmask) == length(stat), stat=$stat"
      end
      # pick sibling with most overlap in down msgs from parent
       # list of similar length siblings
      candidates = dwninters .== maximum(dwninters)
      if candidates[sibidx]
        # must also pick minimized intersect with other remaing siblings
        maxcan = collect(1:len)[candidates]
        with_logger(logger) do
          @info "getSiblingsDelayOrder -- candidates=$candidates, maxcan=$maxcan, rows=$rows"
        end
        if rows[sibidx] == minimum(rows[maxcan])
          with_logger(logger) do
            @info "getSiblingsDelayOrder -- FORCE DOWN INIT SOLVE ON THIS CLIQUE: $(cliq.index), $(cliq.attributes["label"])"
          end
          return false
        end
      end
      with_logger(logger) do
        @info "getSiblingsDelayOrder -- not a max and should block"
      end
      return true
    end

    # still busy solving on branches, so potential to delay
    for i in 1:len
      stillbusymask[i] = maskcol[i] && !(stat[i] in solvedstats)
    end
    with_logger(logger) do
        @info "getSiblingsDelayOrder -- busy solving:"
        @info "maskcol=$maskcol"
        @info "stillbusy=$stillbusymask"
    end

    # Too blunt -- should already have returned false by this point perhaps
    if sum(stillbusymask) > 0
      # yes something to delay about
      with_logger(logger) do
        @info "getSiblingsDelayOrder -- yes delay,"
        @info "stat=$stat"
        @info "symm=$symm"
      end
      return true
    end
  end

  with_logger(logger) do
    @info "getSiblingsDelayOrder -- default will not delay"
  end
  flush(logger.stream)
  # carry over default from partial init process
  return false
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
function getCliqSiblingsPartialNeeds(tree::BayesTree, cliq::Graphs.ExVertex, prnt, dwinmsgs::Dict; logger=ConsoleLogger())
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

  with_logger(logger) do
    @info "getCliqSiblingsPartialNeeds -- CHECK PARTIAL"
  end
  # identify which cliques might have useful information
  localsep = getCliqSeparatorVarIds(cliq)
  seps = Dict{Int, Vector{Symbol}}()
  for si in sibs
    # @show si.attributes["label"]
    mighthave = intersect(getCliqSeparatorVarIds(si), localsep)
    if length(mighthave) > 0
      seps[si.index] = mighthave
      if getCliqStatus(si) in [:initialized; :null; :needdownmsg]
        # partials treated special -- this is slightly hacky
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

  # default
  return false
end

"""
    $SIGNATURES

Bump a clique state machine solver condition in case a task might be waiting on it.
"""
notifyCSMCondition(tree::BayesTree, frsym::Symbol) = notify(getSolveCondition(whichCliq(tree, frsym)))


"""
    $SIGNATURES

Store/cache a clique's solvable dimensions.
"""
function updateCliqSolvableDims!(cliq::Graphs.ExVertex,
                                 sdims::Dict{Symbol, Float64},
                                 logger=ConsoleLogger() )::Nothing
  #
  cliqd = getData(cliq)
  if isready(cliqd.solvableDims)
    take!(cliqd.solvableDims)
    with_logger(logger) do
      @info "cliq $(cliq.index), updateCliqSolvableDims! -- cleared solvableDims"
    end
  end
  put!(cliqd.solvableDims, sdims)
  with_logger(logger) do
      @info "cliq $(cliq.index), updateCliqSolvableDims! -- updated"
  end
  nothing
end

"""
    $SIGNATURES

Retrieve a clique's cached solvable dimensions (since last update).
"""
function fetchCliqSolvableDims(cliq::Graphs.ExVertex)::Dict{Symbol,Float64}
  cliqd = getData(cliq)
  if isready(cliqd.solvableDims)
    return cliqd.solvableDims.data[1]
  end
  return fetch(cliqd.solvableDims)
  # if isready(cliqd.solvableDims)
  # end
end

"""
    $SIGNATURES

Return true there is no other sibling that will make progress.

Notes
- Relies on sibling priority order with only one "currently best" option that will force progress in global upward inference.
- Return false if one of the siblings is still busy
"""
function areSiblingsRemaingNeedDownOnly(tree::BayesTree,
                                        cliq::Graphs.ExVertex  )::Bool
  #
  stillbusylist = [:null; :initialized;]
  prnt = getParent(tree, cliq)
  if length(prnt) > 0
    for si in getChildren(tree, prnt[1])
      # are any of the other siblings still busy?
      if si.index != cliq.index && getCliqStatus(si) in stillbusylist
        return false
      end
    end
  end

  # nope, everybody is waiting for something to change -- proceed with forcing a cliq solve
  return true
end


"""
    $SIGNATURES

Get the main factor graph stored in a history object from a particular step in CSM.

Related

getCliqSubgraphFromHistory, printCliqHistorySummary, printCliqSummary
"""
getGraphFromHistory(hist::Vector{<:Tuple}, step::Int) = hist[step][4].dfg

"""
    $SIGNATURES

Get the cliq sub graph fragment stored in a history object from a particular step in CSM.

Related

getGraphFromHistory, printCliqHistorySummary, printCliqSummary
"""
getCliqSubgraphFromHistory(hist::Vector{<:Tuple}, step::Int) = hist[step][4].cliqSubFg
getCliqSubgraphFromHistory(tree::BayesTree, frnt::Symbol, step::Int) = getCliqSubgraphFromHistory(getCliqSolveHistory(tree, frnt), step)


function printCliqSummary(dfg::G,
                          tree::BayesTree,
                          frs::Symbol,
                          logger=ConsoleLogger() ) where G <: AbstractDFG
  #
  printCliqSummary(dfg, getCliq(tree, frs), logger)
end



function updateSubFgFromDownMsgs!(sfg::G,
                                  dwnmsgs::Dict,
                                  seps::Vector{Symbol} ) where G <: AbstractDFG
  #
  # sanity check basic Bayes (Junction) tree property
  # length(setdiff(keys(dwnmsgs), seps)) == 0 ? nothing : error("updateSubFgFromDownMsgs! -- separators and dwnmsgs not consistent")

  # update specific variables in sfg from msgs
  for (key,beldim) in dwnmsgs
    if key in seps
      setValKDE!(sfg, key, beldim[1], false, beldim[2])
    end
  end

  nothing
end




"""
    $SIGNATURES

Find all factors that go `from` variable to any other complete variable set within `between`.

Notes
- Developed for downsolve in CSM, expanding the cliqSubFg to include all frontal factors.
"""
function findFactorsBetweenFrom(dfg::G,
                                between::Vector{Symbol},
                                from::Symbol ) where {G <: AbstractDFG}
  # get all associated factors
  allfcts = ls(dfg, from)

  # remove candidates with neighbors outside between with mask
  mask = ones(Bool, length(allfcts))
  i = 0
  for fct in allfcts
    i += 1
    # check if immediate neighbors are all in the `between` list
    immnei = ls(dfg, fct)
    if length(immnei) != length(intersect(immnei, between))
      mask[i] = false
    end
  end

  # return only masked factors
  return allfcts[mask]
end

"""
    $SIGNATURES

Special function to add a few variables and factors to the clique sub graph required for downward solve in CSM.

Dev Notes
- There is still some disparity on whether up and down solves of tree should use exactly the same subgraph...  'between for up and frontal connected for down'
"""
function addDownVariableFactors!(dfg::G1,
                                 subfg::G2,
                                 cliq::Graphs.ExVertex,
                                 logger=ConsoleLogger()  ) where {G1 <: AbstractDFG, G2 <: InMemoryDFGTypes}
  #
  # determine which variables and factors needs to be added
  currsyms = ls(subfg)
  allclsyms = getCliqVarsWithFrontalNeighbors(dfg, cliq)
  newsyms = setdiff(allclsyms, currsyms)
  with_logger(logger) do
    @info "addDownVariableFactors!, cliq=$(cliq.index), newsyms=$newsyms"
  end
  frtls = getCliqFrontalVarIds(cliq)
  with_logger(logger) do
    @info "addDownVariableFactors!, cliq=$(cliq.index), frtls=$frtls"
  end
  allnewfcts = union(map(x->findFactorsBetweenFrom(dfg,union(currsyms, newsyms),x), frtls)...)
  newfcts = setdiff(allnewfcts, lsf(subfg))
  with_logger(logger) do
    @info "addDownVariableFactors!, cliq=$(cliq.index), newfcts=$newfcts, allnewfcts=$allnewfcts"
  end

  # add the variables
  DFG.getSubgraph(dfg, newsyms, false, subfg)
  # add the factors
  DFG.getSubgraph(dfg, newfcts, false, subfg)

  return newsyms, newfcts
end


"""
    $SIGNATURES

Basic wrapper to take local product and then set the value of `sym` in `dfg`.
"""
function localProductAndUpdate!(dfg::G,
                                sym::Symbol,
                                setkde::Bool=true,
                                logger=ConsoleLogger() )::Tuple{BallTreeDensity, Float64, Vector{Symbol}} where {G <: AbstractDFG}
  #
  pp, dens, parts, lbl, infdim = localProduct(dfg, sym, logger=logger)
  setkde ? setValKDE!(dfg, sym, pp, false, infdim) : nothing

  return pp, infdim, lbl
end

"""
    $SIGNATURES

Calculate new and then set the down messages for a clique in Bayes (Junction) tree.
"""
function getSetDownMessagesComplete!(subfg::G,
                                     cliq::Graphs.ExVertex,
                                     prntDwnMsgs::TempBeliefMsg,
                                     logger=ConsoleLogger()  )::Nothing where G <: AbstractDFG
  #
  allvars = getCliqVarIdsAll(cliq)
  allprntkeys = collect(keys(prntDwnMsgs))
  passkeys = intersect(allvars, setdiff(allprntkeys,ls(subfg)))
  remainkeys = setdiff(allvars, passkeys)
  newDwnMsgs = TempBeliefMsg()

  # some msgs are just pass through from parent
  for pk in passkeys
    newDwnMsgs[pk] = prntDwnMsgs[pk]
  end

  # other messages must be extracted from subfg
  for mk in remainkeys
    newDwnMsgs[mk] = (getKDE(subfg, mk), getVariableInferredDim(subfg,mk))
  end

  # set the downward keys
  with_logger(logger) do
    @info "cliq $(cliq.index), getSetDownMessagesComplete!, allkeys=$(allvars), passkeys=$(passkeys)"
  end
  setDwnMsg!(cliq, newDwnMsgs)

  return nothing
end


"""
    $SIGNATURES

Determine which variables to iterate or compute directly for downward tree pass of inference.

Related Functions from Upward Inference

directPriorMsgIDs, directFrtlMsgIDs, directAssignmentIDs, mcmcIterationIDs
"""
function determineCliqVariableDownSequence(subfg::G, cliq::Graphs.ExVertex) where G <: AbstractDFG
  frtl = getCliqFrontalVarIds(cliq)

  adj = DFG.getAdjacencyMatrix(subfg)
  mask = map(x->(x in frtl), adj[1,:])
  subAdj = adj[2:end,mask] .!= nothing
  newFrtlOrder = Symbol.(adj[1,mask])
  crossCheck = sum(Int.(subAdj), dims=2) .> 1
  iterVars = Symbol[]
  for i in 1:length(crossCheck)
    # must add associated variables to iterVars
    if crossCheck[i]
      # find which variables are associated
      varSym = newFrtlOrder[subAdj[i,:]]
      union!(iterVars, varSym)
    end
  end

  # return iteration list ordered by frtl
  return intersect(frtl, iterVars)
end


"""
    $SIGNATURES

Perform downward direction solves on a sub graph fragment.
Calculates belief on each of the frontal variables and iterate if required.

Notes
- uses all factors connected to the frontal variables.
- assumes `subfg` was properly prepared before calling.
- has multi-process option.

Dev Notes
- TODO incorporate variation possible due to cross frontal factors.
- cleanup and updates required, and @spawn jl 1.3
"""
function solveCliqDownFrontalProducts!(subfg::G,
                                       cliq::Graphs.ExVertex,
                                       opts::SolverParams,
                                       logger=ConsoleLogger();
                                       MCIters::Int=3 )::Nothing where G <: AbstractDFG
  #
  # get frontal variables for this clique
  frsyms = getCliqFrontalVarIds(cliq)

  # determine if cliq has cross frontal factors
  # iterdwn, directdwns, passmsgs?
  iterFrtls = determineCliqVariableDownSequence(subfg,cliq)

  # direct frontals
  directs = setdiff(frsyms, iterFrtls)

  # ignore limited fixed lag variables
  fixd = map(x->opts.limitfixeddown && isMarginalized(subfg,x), frsyms)
  skip = frsyms[fixd]
  iterFrtls = setdiff(iterFrtls, skip)
  directs = setdiff(directs, skip)
  with_logger(logger) do
    @info "cliq $(cliq.index), doCliqDownSolve_StateMachine, skipping marginalized keys=$(skip)"
  end


  # use new localproduct approach
  if opts.multiproc
    downresult = Dict{Symbol, Tuple{BallTreeDensity, Float64, Vector{Symbol}}}()
    @sync for i in 1:length(directs)
      @async begin
        downresult[directs[i]] = remotecall_fetch(localProductAndUpdate!, upp2(), subfg, directs[i], false)
      end
    end
    with_logger(logger) do
      @info "cliq $(cliq.index), doCliqDownSolve_StateMachine, multiproc keys=$(keys(downresult))"
    end
    for fr in directs
        with_logger(logger) do
            @info "cliq $(cliq.index), doCliqDownSolve_StateMachine, key=$(fr), infdim=$(downresult[fr][2]), lbls=$(downresult[fr][3])"
        end
      setValKDE!(subfg, fr, downresult[fr][1], false, downresult[fr][2])
    end
    for mc in 1:MCIters, fr in iterFrtls
      result = remotecall_fetch(localProductAndUpdate!, upp2(), subfg, fr, false)
      setValKDE!(subfg, fr, result[1], false, result[2])
      with_logger(logger) do
          @info "cliq $(cliq.index), doCliqDownSolve_StateMachine, iter key=$(fr), infdim=$(result[2]), lbls=$(result[3])"
      end
    end
  else
    # do directs first
    for fr in directs
      localProductAndUpdate!(subfg, fr, true, logger)
    end
    #do iters next
    for mc in 1:MCIters, fr in iterFrtls
      localProductAndUpdate!(subfg, fr, true, logger)
    end
  end

  return nothing
end




#
