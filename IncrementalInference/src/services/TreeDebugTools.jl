
#

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
function getCliqSubgraphFromHistory(
  tree::AbstractBayesTree,
  hists::Dict{Symbol, Tuple},
  frnt::Symbol,
  step::Int,
)
  return getCliqSubgraphFromHistory(hists[frnt], step)
end

function printCliqSummary(
  dfg::G,
  tree::AbstractBayesTree,
  frs::Symbol,
  logger = ConsoleLogger(),
) where {G <: AbstractDFG}
  #
  return printCliqSummary(dfg, getClique(tree, frs), logger)
end

"""
    $(SIGNATURES)

Calculate a fresh (single step) approximation to the variable `sym` in clique `cliq` as though during the upward message passing.  The full inference algorithm may repeatedly calculate successive apprimxations to the variables based on the structure of the clique, factors, and incoming messages.
Which clique to be used is defined by frontal variable symbols (`cliq` in this case) -- see `getClique(...)` for more details.  The `sym` symbol indicates which symbol of this clique to be calculated.  **Note** that the `sym` variable must appear in the clique where `cliq` is a frontal variable.
"""
function treeProductUp(
  fg::AbstractDFG,
  tree::AbstractBayesTree,
  cliq::Symbol,
  sym::Symbol;
  N::Int = 100,
  dbg::Bool = false,
)
  #
  cliq = getClique(tree, cliq)
  cliqdata = getCliqueData(cliq)

  # perform the actual computation
  pGM, fulldim = predictbelief(fg, sym, :; N = N, dbg = dbg)

  return pGM, nothing
end

"""
    $(SIGNATURES)

Calculate a fresh---single step---approximation to the variable `sym` in clique `cliq` as though during the downward message passing.  The full inference algorithm may repeatedly calculate successive apprimxations to the variable based on the structure of variables, factors, and incoming messages to this clique.
Which clique to be used is defined by frontal variable symbols (`cliq` in this case) -- see `getClique(...)` for more details.  The `sym` symbol indicates which symbol of this clique to be calculated.  **Note** that the `sym` variable must appear in the clique where `cliq` is a frontal variable.
"""
function treeProductDwn(
  fg::G,
  tree::AbstractBayesTree,
  cliq::Symbol,
  sym::Symbol;
  N::Int = 100,
  dbg::Bool = false,
) where {G <: AbstractDFG}
  #
  @warn "treeProductDwn might not be working properly at this time. (post DFG v0.6 upgrade maintenance required)"
  cliq = getClique(tree, cliq)
  cliqdata = getCliqueData(cliq)

  # get the local variable id::Int identifier
  # vertid = fg.labelDict[sym]

  # get all the incoming (upward) messages from the tree cliques
  # convert incoming messages to Int indexed format (semi-legacy format)
  cl = parentCliq(tree, cliq)
  msgdict = getDwnMsgs(cl[1])
  dict = Dict{Int, TreeBelief}()
  for (dsy, btd) in msgdict
    dict[fg.IDs[dsy]] =
      TreeBelief(btd.val, btd.bw, btd.infoPerCoord, getStateKind(getVariable(fg, sym)))
  end
  dwnmsgssym = LikelihoodMessage[LikelihoodMessage(dict);]

  # perform the actual computation
  pGM, fulldim = predictbelief(fg, sym, :; N = N, dbg = dbg)

  return pGM, nothing, sym, dwnmsgssym
end

"""
    $SIGNATURES

Print one specific line of a clique state machine history.

Related:

[`printCliqHistorySequential`](@ref), [`printCliqHistorySummary`](@ref)
"""
function printHistoryLine(
  fid,
  hi::CSMHistoryTuple,
  cliqid::AbstractString = "",
  seq::Int = 0,
)
  #

  # global sequence number
  first = clampBufferString("$seq", 5)

  # 5.13
  first *= clampBufferString("$cliqid.$(string(hi[2]))", 6)

  # time
  first *= clampBufferString(string(split(string(hi[1]), 'T')[end]), 14)

  # next function
  nextfn = split(split(string(hi[3]), '.')[end], '_')[1]
  first *= clampBufferString(nextfn * "                             ", 20, 18)

  first *= " | "
  first *= clampBufferString(string(getCliqueStatus(hi[4].cliq)) * "             ", 9, 7)

  # parent status
  first *= " P "
  downRxMessage = getMessageBuffer(hi[4].cliq).downRx
  toadd = if !isnothing(downRxMessage)
    #TODO maybe don't use tree here
    "$(getParent(hi[4].tree, hi[4].cliq)[1].id):$(downRxMessage.status)"
  else
    "  ----"
  end
  first *= clampBufferString(toadd * "          ", 9, 7)

  # children status
  first = first * "C "

  upRxMessages = getMessageBuffer(hi[4].cliq).upRx
  # all_child_status = map((k,msg) -> (k,msg.status), pairs(upRxMessages))
  if length(upRxMessages) > 0
    for (k, msg) in upRxMessages
      toadd = string(k) * ":" * string(msg.status) * " "
      first *= clampBufferString(toadd * "                  ", 8, 7)
    end
  else
    first *= clampBufferString("  ----", 8)
  end

  # sibling status # TODO JT removed but kept for future if needed
  # first *= "|S| "
  # if 0 < length(hi[4].parentCliq)
  #   frt = (hi[4].parentCliq[1] |> getFrontals)[1]
  #   childs = getChildren(hi[4].tree, frt)
  #   # remove current clique to leave only siblings
  #   filter!(x->x.index!=hi[4].cliq.id.value, childs)
  #   for ch in childs
  #     first = first*"$(ch.index)"*string(getCliqueStatus(ch))*" "
  #   end
  # end

  return println(fid, first)
end

"""
    $SIGNATURES

Print a short summary of state machine history for a clique solve.

Related:

getTreeAllFrontalSyms, animateCliqStateMachines, printHistoryLine, printCliqHistorySequential
"""
function printCliqHistorySummary(
  fid,
  hist::Vector{CSMHistoryTuple},
  cliqid::AbstractString = "",
)
  if length(hist) == 0
    @warn "printCliqHistorySummary -- No CSM history found."
  end
  for hi in hist
    printHistoryLine(fid, hi, cliqid)
  end
  return nothing
end

function printCliqHistorySummary(hist::Vector{CSMHistoryTuple}, cliqid::AbstractString = "")
  #
  return printCliqHistorySummary(stdout, hist, cliqid)
end

function printCliqHistorySummary(
  hists::Dict{Int, Vector{CSMHistoryTuple}},
  tree::AbstractBayesTree,
  sym::Symbol,
)
  #
  hist = hists[getClique(tree, sym).id]
  return printCliqHistorySummary(stdout, hist, string(getClique(tree, sym).id))
end

# TODO maybe Base. already has something like this Union{UnitRange, AbstractVector, etc.}
const CSMRangesT{T} = Union{T, UnitRange{T}, <:AbstractVector{T}}
const CSMRanges = CSMRangesT{Int}
# old
# const CSMTupleRangesT{T} = Union{Tuple{T,T},Tuple{T,UnitRange{T}},Tuple{T,AbstractVector{T}},Tuple{UnitRange{T},T},Tuple{UnitRange{T},UnitRange{T}},Tuple{AbstractVector{T},T},Tuple{AbstractVector{T},AbstractVector{T}},Tuple{AbstractVector{T},UnitRange{T}} }

"""
    $SIGNATURES

Print a sequential summary lines of clique state machine histories in hists::Dict.

Notes
- Slices are allowed, see examples.

Example
```julia
printCSMHistorySequential(hists)
printCSMHistorySequential(hists, 2=>46)
printCSMHistorySequential(hists, 1=>11:15)
printCSMHistorySequential(hists, [1,4,6]=>11:15)
printCSMHistorySequential(hists, [2=>45:52; 1=>10:15])
```

DevNotes
- TODO perhaps move some of this functionality upstream to FSM
- TODO upgrade to default `whichstep = :=>:` -- i.e. 
  - add dispatch for `(:) |> typeof == Colon`, 
  - `5:6=>:`.
- TODO also add a elements between `Tuple{<:CSMRanges, Pair{<:CSMRanges,<:CSMRanges}}` option
  - ([1;3], 1=>5:7) which will print all steps from CSM 1 and 3, which occur between 1=>5 and 1=>7.
- TODO maybe also `Dict(5=>5:8, 8=>20:25)`, or `Dict(2:5=>[1;3], 10=>1:5)`.

Related:

printHistoryLine, printCliqHistory
"""
function printCSMHistorySequential(
  hists::Dict{Int, Vector{CSMHistoryTuple}},
  whichsteps::Union{Nothing, Vector{<:Pair{<:CSMRanges, <:CSMRanges}}} = nothing,
  fid = stdout,
)
  #
  # vectorize all histories in single Array
  allhists = Vector{CSMHistoryTuple}()
  alltimes = Vector{DateTime}()
  allcliqids = Vector{Int}()
  for (cid, hist) in hists, hi in hist
    push!(allhists, hi)
    push!(alltimes, hi[1])
    push!(allcliqids, cid)
  end

  # sort array by timestamp element
  pm = sortperm(alltimes)
  allhists_ = allhists[pm]
  alltimes_ = alltimes[pm]
  allcliqids_ = allcliqids[pm]

  # print each line of the sorted array with correct cliqid marker
  for idx = 1:length(alltimes)
    hiln = allhists_[idx]
    # show only one line if whichstep is not nothing
    inSliceList = whichsteps === nothing
    if !inSliceList
      for whichstep in whichsteps
        inSliceList && break
        inSliceList =
          inSliceList || (allcliqids_[idx] in whichstep[1] && hiln[2] in whichstep[2])
      end
    end
    if inSliceList
      printHistoryLine(fid, hiln, string(allcliqids_[idx]), idx)
    end
  end
  return nothing
end

function printCSMHistorySequential(
  hists::Dict{Int, Vector{CSMHistoryTuple}},
  whichstep::Pair{<:CSMRanges, <:CSMRanges},
  fid = stdout,
)
  #
  return printCSMHistorySequential(hists, [whichstep;], fid)
end

function printCSMHistorySequential(
  hists::Dict{Int, Vector{CSMHistoryTuple}},
  fid::AbstractString,
)
  #
  @info "printCliqHistorySequential -- assuming file request, writing history to $fid"
  file = open(fid, "w")
  printCSMHistorySequential(hists, nothing, file)
  close(file)
  return nothing
end

"""
    $SIGNATURES

Print one line of lanes summarizing all clique state machine histories.

Notes
- hiVec is vector of all cliques (i.e. lanes) to print as one LINE into `fid` 
  - contains `::Tuple{Int,..}` with global counter (not the default CSM counter)
- Vector of `CSMHistoryTuple`

Related:

printCliqHistoryLogical, printCliqHistoryLine
"""
function printHistoryLane(
  fid,
  linecounter::Union{Int, String},
  hiVec::Vector{<:Union{NamedTuple, Tuple}},
  seqLookup::NothingUnion{Dict{Pair{Int, Int}, Int}} = nothing,
)
  #

  ## build a string
  line = clampBufferString("$linecounter", 4)
  for counter = 1:length(hiVec)
    # lane marker
    line *= "| "
    if !isassigned(hiVec, counter)
      line *= clampBufferString("", 19)
      continue
    end
    hi = hiVec[counter]
    # global counter
    useCount = seqLookup !== nothing ? seqLookup[(hi[4].cliq.id.value => hi[2])] : hi[2]
    line *= clampBufferString("$(useCount)", 4)
    # next function
    nextfn = split(string(hi[3]), '.')[end]
    line *= clampBufferString(nextfn, 10, 9)
    # clique status
    st = hi[4] isa String ? hi[4] : string(getCliqueStatus(hi[4].cliq))
    line *= clampBufferString(st, 5, 4)
  end
  ## print the string
  return println(fid, line)
end

"""
    $SIGNATURES

Print history in swimming lanes, side by side with global sequence counter.

Examples

```julia
printCSMHistoryLogical(hist)
printCSMHistoryLogical(hists, order=[4;3], printLines=2:5)

# or to a IOStream, file, network, etc
fid = open(joinLogPath(fg, "CSMHistCustom.txt"),"w")
printCSMHistoryLogical(hist, fid)
close(fid)
```

DevNotes
- `order` should be flexible like `Sequential` and `<:CSMRanges`.
"""
function printCSMHistoryLogical(
  hists::Dict{Int, Vector{CSMHistoryTuple}},
  fid = stdout;
  order::AbstractVector{Int} = sort(collect(keys(hists))),
  printLines = 1:99999999,
)
  #

  # vectorize all histories in single Array
  allhists = Vector{CSMHistoryTuple}()
  alltimes = Vector{DateTime}()
  allcliqids = Vector{Int}()
  # "lanes" (i.e. individual cliques)
  numLanes = length(order)
  # "lines" (i.e. CSM steps)
  maxLines = [0; 0]
  for (cid, hist) in hists
    # find max number of lines to print later
    maxLines[2] = length(hist)
    maxLines[1] = maximum(maxLines)
    for hi in hist
      push!(allhists, hi)
      push!(alltimes, hi[1])
      push!(allcliqids, cid)
    end
  end
  maxLines[1] = minimum([maxLines[1]; printLines[end]])

  # sort array by timestamp element
  pm = sortperm(alltimes)
  allhists_ = allhists[pm]
  alltimes_ = alltimes[pm]
  allcliqids_ = allcliqids[pm]

  # first get the global sequence (invert order dict as bridge table)
  seqLookup = Dict{Pair{Int, Int}, Int}()
  for idx = 1:length(alltimes)
    hiln = allhists_[idx]
    seqLookup[(hiln[4].cliq.id.value => hiln[2])] = idx
  end

  # print the column titles
  titles = Vector{Tuple{String, Int, String, String}}()
  for ord in order
    csym = 0 < length(hists[ord]) ? getFrontals(hists[ord][1][4].cliq)[1] |> string : ""
    csym = clampBufferString("$csym", 9)
    push!(titles, ("", ord, csym, clampBufferString("", 10)))
  end
  printHistoryLane(fid, "", titles)
  print(fid, "----")
  for i = 1:numLanes
    print(fid, "+--------------------")
  end
  println(fid, "")

  glbSeqCt = 0 # Ref{Int}(0)
  ## repeat for the maximum number of "lines" (i.e. CSM steps)
  for idx = printLines[1]:maxLines[1]

    ## build each line as vector of "lanes" (i.e. individual cliques)
    allLanes = Vector{CSMHistoryTuple}(undef, numLanes)

    laIdx = 0
    for laId in order
      laIdx += 1
      # if history data exists for this line (idx) and lane (laId), then build a new lane Tuple
      if idx <= length(hists[laId])
        # FIXME, swat first element with global counter (not local as stored in hists)
        # @show hists[laId][idx][3]
        allLanes[laIdx] = hists[laId][idx]
      end
    end

    #$cliqid.$(string(hi[2]))
    # glbSeqCt += 1
    printHistoryLane(fid, idx, allLanes, seqLookup)
  end
end
# print each line of the sorted array with correct cliqid marker

"""
  $SIGNATURES

Repeat a solver state machine step -- useful for debugging. 

Notes
- use in combination with `solveTree!(fg, recordcliqs=[:0; :x7; ...])` -- i.e. by clique frontals as identifier
  - to record everything, one can do: `recordcliqs=ls(fg)`.
- `duplicate` avoids changing history or prime data in `hists`.
- Replaces old API `sandboxCliqResolveStep`
- Consider using this in combination with tools like [Revise.jl](https://github.com/timholy/Revise.jl)
  - On by default in VSCode.
- Internally sets `csmc.enableLogging=false`

Example

```julia
using IncrementalInference

# generate a factor graph
fg = generateGraph_Kaess()

# solve and record everything
smtasks = Task[]
tree, _, = solveTree!(fg, smtasks=smtasks, recordcliqs=ls(fg));
# draw Bayes tree with graphviz and xdot installed
drawTree(tree, show=true)

# fetch histories
hists = fetchCliqHistoryAll!(smtasks);

# check a new csmc before step 2
csmc_ = repeatCSMStep!(hists, 1, 1)

# For use with VSCode debugging
@enter repeatCSMStep!(hists, 1, 1)

# or perhaps test a longer chain of changes
hists_ = deepcopy(hists)
repeatCSMStep!(hists_, 1, 4, duplicate=false)
repeatCSMStep!(hists_, 1, 5, duplicate=false)
repeatCSMStep!(hists_, 1, 6, duplicate=false)
```

DevNotes
- TODO consolidate upstream with `FSM.sandboxStateMachineStep`

Related

[`solveTree!`](@ref), [`solveCliqUp!`](@ref), [`fetchCliqHistoryAll`](@ref), [`printCSMHistoryLogical`](@ref), [`printCSMHistorySequential`](@ref), cliqHistFilterTransitions
"""
function repeatCSMStep!(
  hist::AbstractVector{<:CSMHistoryTuple},
  step::Int;
  duplicate::Bool = true,
  enableLogging::Bool = false,
)
  #

  # the function at steo 
  fnc_ = hist[step].f
  # the data before step
  csmc_ = (duplicate ? x -> deepcopy(x) : x -> x)(hist[step].csmc)
  csmc_.enableLogging = enableLogging
  csmc_.logger = enableLogging ? SimpleLogger() : SimpleLogger(Base.devnull)

  # run the step
  newfnc_ = fnc_(csmc_)
  return newfnc_, csmc_
end

function repeatCSMStep!(
  hists::Dict{Int, <:AbstractVector{CSMHistoryTuple}},
  csmid::Int,
  step::Int;
  duplicate::Bool = true,
  enableLogging::Bool = false,
)
  #

  return repeatCSMStep!(
    hists[csmid],
    step;
    duplicate = duplicate,
    enableLogging = enableLogging,
  )
end

"""
    $SIGNATURES
Reattach a CSM's data container after the deepcopy used from recordcliq.

MIGHT BE OBSOLETE
"""
function attachCSM!(
  csmc::CliqStateMachineContainer,
  dfg::AbstractDFG,
  tree::MetaBayesTree;
  logger = SimpleLogger(),
)
  #
  # csmc = csmc__

  @error("attachCSM! has been updated without testing and might not work as you intended.")

  csmc.dfg = dfg
  csmc.tree = tree
  csmc.logger = logger # TODO option to reopen and append to previous logger file

  @info "attaching csmc and dropping any contents from csmc's previously held (copied) message channels."
  cid = csmc.cliq.id.value
  # pids = csmc.parentCliq .|> x->x.id
  # cids = csmc.childCliqs .|> x->x.id

  csmc.cliq = getClique(tree, cid)
  # csmc.parentCliq = pids .|> x->getindex(tree.cliques, x)
  # csmc.childCliqs = cids .|> x->getindex(tree.cliques, x)

  return csmc
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
function animateCliqStateMachines(
  tree::AbstractBayesTree,
  cliqsyms::Vector{Symbol},
  hists::Dict{Symbol, Tuple};
  frames::Int = 100,
)
  #
  error("`animateCliqStateMachines` is outdated")
  startT = Dates.now()
  stopT = Dates.now()

  # get start and stop times across all cliques
  first = true
  for sym in cliqsyms
    hist = hists[sym] #getCliqSolveHistory(tree, sym)
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
  folders = String[]
  for sym in cliqsyms
    hist = hists[sym] #getCliqSolveHistory(tree, sym)
    # hist = getCliqSolveHistory(tree, sym)
    retval = animateStateMachineHistoryByTime(
      hist;
      frames = frames,
      folder = "caesar/animatecsm/cliq$sym",
      title = "$sym",
      startT = startT,
      stopT = stopT,
      rmfirst = false,
    )
    push!(folders, "cliq$sym")
  end

  return folders
end

"""
    $SIGNATURES

Return state machine transition steps from history such that the `nextfnc::Function`.

Related:

printCliqHistorySummary, filterHistAllToArray, sandboxCliqResolveStep
"""
function cliqHistFilterTransitions(hist::Vector{CSMHistoryTuple}, nextfnc::Function)
  ret = Vector{CSMHistoryTuple}()
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

printCliqHistorySummary, cliqHistFilterTransitions, sandboxCliqResolveStep
"""
function filterHistAllToArray(
  tree::AbstractBayesTree,
  hists::Dict{Symbol, Tuple},
  frontals::Vector{Symbol},
  nextfnc::Function,
)
  error("filterHistAllToArray needs to be updated for new CSM")
  ret = Vector{CSMHistoryTuple}()
  for sym in frontals
    hist = hists[sym] # getCliqSolveHistory(tree, sym)
    fih = cliqHistFilterTransitions(hist, nextfnc)
    for fi in fih
      push!(ret, fi)
    end
  end
  return ret
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
tree = solveTree!(fg, recordcliqs=fsy)

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
function animateCSM(
  tree::AbstractBayesTree,
  autohist::Dict{Int, T};
  frames::Int = 100,
  interval::Int = 2,
  dpi::Int = 100,
  rmfirst::Bool = true,
  folderpath::AbstractString = "/tmp/caesar/csmCompound/",
  fsmColors::Dict{Symbol, String} = Dict{Symbol, String}(),
  defaultColor::AbstractString = "red",
) where {T <: AbstractVector}
  #

  easyNames = Dict{Symbol, Int}()
  hists =
    Dict{Symbol, Vector{Tuple{DateTime, Int64, Function, CliqStateMachineContainer}}}()
  for (id, hist) in autohist
    frtl = getFrontals(getClique(tree, id))
    hists[frtl[1]] = Vector{Tuple{DateTime, Int64, Function, CliqStateMachineContainer}}()
    for hi in hist
      push!(hists[frtl[1]], (hi.timestamp, hi.id, hi.f, hi.csmc)) # Tuple.(hist)
    end
    easyNames[frtl[1]] = id
  end

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
    @warn "Removing $folderpath in preparation for new frames."
    Base.rm("$folderpath"; recursive = true, force = true)
  end

  function csmTreeAni(hl::Tuple, frame::Int, folderpath::AbstractString)
    drawTree(
      hl[4].tree;
      show = false,
      filepath = joinpath(folderpath, "tree_$frame.png"),
      dpi = dpi,
    )
    return nothing
  end

  function autocolor_cb(hi::Tuple, csym::Symbol, aniT::DateTime)
    retc = getCliqueDrawColor(hi[4].cliq)
    return (retc === nothing ? "gray" : retc)
  end

  # animateStateMachineHistoryByTimeCompound(hists, startT, stopT, folder="caesar/csmCompound", frames=frames)
  return FSM.animateStateMachineHistoryIntervalCompound(
    hists;
    easyNames = easyNames,
    folderpath = folderpath,
    interval = interval,
    dpi = dpi,
    draw_more_cb = csmTreeAni,
    fsmColors = fsmColors,
    defaultColor = defaultColor,
    autocolor_cb = autocolor_cb,
  )
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
function makeCsmMovie(
  fg::AbstractDFG,
  tree::AbstractBayesTree,
  cliqs = ls(fg);
  assignhist = nothing,
  show::Bool = true,
  filename::AbstractString = "/tmp/caesar/csmCompound/out.ogv",
  frames::Int = 1000,
)
  #
  if assignhist != nothing
    assignTreeHistory!(tree, assignhist)
  end
  csmAnimate(fg, tree, cliqs; frames = frames)
  # Base.rm("/tmp/caesar/csmCompound/out.ogv")
  run(
    `ffmpeg -r 10 -i /tmp/caesar/csmCompound/csm_%d.png -c:v libtheora -vf fps=25 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -q 10 $filename`,
  )
  if show
    @async run(`totem $filename`)
  end
  return filename
end

#
