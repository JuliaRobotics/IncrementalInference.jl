
### DONT DELETE YET -- see more likely list below


# function getShortestPathNeighbors(fgl::FactorGraph;
#     from::TreeClique=nothing,
#     to::TreeClique=nothing,
#     neighbors::Int=0 )
#
#   edgelist = shortest_path(fgl.g, ones(num_edges(fgl.g)), from, to)
#   vertdict = Dict{Int,TreeClique}()
#   edgedict = edgelist2edgedict(edgelist)
#   expandVertexList!(fgl, edgedict, vertdict) # grow verts
#   for i in 1:neighbors
#     expandEdgeListNeigh!(fgl, vertdict, edgedict) # grow edges
#     expandVertexList!(fgl, edgedict, vertdict) # grow verts
#   end
#   return vertdict
# end

# function subgraphShortestPath(fgl::FactorGraph;
#                               from::TreeClique=nothing,
#                               to::TreeClique=nothing,
#                               neighbors::Int=0  )
#   #
#   vertdict = getShortestPathNeighbors(fgl, from=from, to=to, neighbors=neighbors)
#   return genSubgraph(fgl, vertdict)
# end



"""
    $SIGNATURES

Draw and show the factor graph `<:AbstractDFG` via system graphviz and pdf app.

Notes
- Should not be calling outside programs.
- Need long term solution
- DFG's `toDotFile` a better solution -- view with `xdot` application.
- also try `engine={"sfdp","fdp","dot","twopi","circo","neato"}`

Future:
- Might be kept with different call strategy since this function so VERY useful!
- Major issue that this function calls an external program such as "evince", which should be
   under user control only.
- Maybe solution is
- ```toDot(fg,file=...); @async run(`xdot file.dot`)```, or
  - ```toDot(fg,file=...); exportPdf(...); @async run(`evince ...pdf`)```.
"""
function writeGraphPdf(fgl::G;
                       viewerapp::String="evince",
                       filepath::AS="/tmp/fg.pdf",
                       engine::AS="neato", #sfdp
                       show::Bool=true ) where {G <: AbstractDFG, AS <: AbstractString}
  #
  @warn "writeGraphPdf is function changing to drawGraph, see DFG.toDotFile(dfg) as part of the long term solution."

  fgd = fgl
  @info "Writing factor graph file"
  fext = split(filepath, '.')[end]
  fpwoext = filepath[1:(end-length(fext)-1)] # split(filepath, '.')[end-1]
  dotfile = fpwoext*".dot"

  # create the dot file
  DFG.toDotFile(fgl, dotfile)

  try
    run(`$(engine) $(dotfile) -T$(fext) -o $(filepath)`)
    show ? (@async run(`$(viewerapp) $(filepath)`)) : nothing
  catch e
    @warn "not able to show $(filepath) with viewerapp=$(viewerapp). Exception e=$(e)"
  end
  nothing
end





function downGibbsCliqueDensity(inp::ExploreTreeType{T},
                                N::Int=100,
                                dbg::Bool=false,
                                logger=ConsoleLogger()  ) where {T}
  #
  with_logger(logger) do
    @info "=================== Iter Clique $(getLabel(inp.cliq)) ==========================="
  end
  mcmciter = inp.prnt != nothing ? 3 : 0
  return downGibbsCliqueDensity(inp.fg, inp.cliq, inp.sendmsgs, N, mcmciter, dbg)
end

function prepDwnPreOrderStack!(bt::AbstractBayesTree,
                               parentStack::Array{TreeClique,1})
  # dwn message passing function
  nodedata = nothing
  tempStack = Array{TreeClique,1}()
  push!(tempStack, parentStack[1])

  while ( length(tempStack) != 0 || nodedata != nothing )
    if nodedata != nothing
      for child in out_neighbors(nodedata, bt.bt) #nodedata.cliq, nodedata.bt.bt
          ## TODO -- don't yet have the NBPMessage results to propagate belief
          ## should only construct ett during processing
          push!(parentStack, child) #ett
          push!(tempStack, child) #ett
      end
      nodedata = nothing # its over but we still need to complete computation in each of the leaves of the tree
    else
      nodedata = tempStack[1]
      deleteat!(tempStack, 1)
    end
  end
  nothing
end

function findVertsAssocCliq(fgl::FactorGraph, cliq::TreeClique)

  cliqdata = getData(cliq)
  IDS = [cliqdata.frontalIDs; cliqdata.separatorIDs] #inp.cliq.attributes["frontalIDs"]


  @error "findVertsAssocCliq -- not completed yet"
  nothing
end

function partialExploreTreeType(pfg::G, pbt::AbstractBayesTree, cliqCursor::TreeClique, prnt, pmsgs::Array{NBPMessage,1}) where G <: AbstractDFG
    # info("starting pett")
    # TODO -- expand this to grab only partial subsection from the fg and bt data structures


    if length(pmsgs) < 1
      return ExploreTreeType(pfg, pbt, cliqCursor, prnt, NBPMessage[])
    else
      return ExploreTreeType(pfg, pbt, cliqCursor, prnt, pmsgs)
    end
    nothing
end

function dispatchNewDwnProc!(fg::G,
                             bt::AbstractBayesTree,
                             parentStack::Array{TreeClique,1},
                             stkcnt::Int,
                             refdict::Dict{Int,Future};
                             N::Int=100,
                             dbg::Bool=false,
                             drawpdf::Bool=false  ) where G <: AbstractDFG
  #
  cliq = parentStack[stkcnt]
  while !haskey(refdict, cliq.index) # nodedata.cliq
    sleep(0.25)
  end

  rDDT = fetch(refdict[cliq.index]) #nodedata.cliq
  delete!(refdict, cliq.index) # nodedata

  if rDDT != nothing
    updateFGBT!(fg, bt, cliq.index, rDDT, dbg=dbg, fillcolor="lightblue")
    setCliqStatus!(cliq, :downsolved) # should be a notify
    drawpdf ? drawTree(bt) : nothing
  end

  emptr = BayesTree(nothing, 0, Dict{Int,TreeClique}(), Dict{String,Int}());

  for child in out_neighbors(cliq, bt.bt) # nodedata.cliq, nodedata.bt.bt
      haskey(refdict, child.index) ? error("dispatchNewDwnProc! -- why you already have dwnremoteref?") : nothing
      ett = partialExploreTreeType(fg, emptr, child, cliq, [rDDT.dwnMsg]) # bt
      refdict[child.index] = remotecall(downGibbsCliqueDensity, upp2() , ett, N) # Julia 0.5 swapped order
  end
  nothing
end

"""
    $SIGNATURES

Downward message passing on Bayes (Junction) tree.

Notes
- Simultaenously launches as many async dispatches to remote processes as there are cliques in the tree.
"""
function processPreOrderStack!(fg::G,
                               bt::AbstractBayesTree,
                               parentStack::Array{TreeClique,1},
                               refdict::Dict{Int,Future};
                               N::Int=100,
                               dbg::Bool=false,
                               drawpdf::Bool=false ) where G <: AbstractDFG
  #
    # dwn message passing function for iterative tree exploration
    stkcnt = 0

    @sync begin
      sendcnt = 1:length(parentStack) # separate memory for remote calls
      for i in 1:sendcnt[end]
          @async try
            dispatchNewDwnProc!(fg, bt, parentStack, sendcnt[i], refdict, N=N, dbg=dbg, drawpdf=drawpdf) # stkcnt ##pidxI,nodedata
          catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)
          end
      end
    end
    nothing
end

function downMsgPassingIterative!(startett::ExploreTreeType{T};
                                  N::Int=100,
                                  dbg::Bool=false,
                                  drawpdf::Bool=false  ) where {T}
  #
  # this is where we launch the downward iteration process from
  parentStack = Array{TreeClique,1}()
  refdict = Dict{Int,Future}()

  # start at the given clique in the tree -- shouldn't have to be the root.
  pett = partialExploreTreeType(startett.fg, startett.bt, startett.cliq,
                                        startett.prnt, startett.sendmsgs)
  refdict[startett.cliq.index] = remotecall(downGibbsCliqueDensity, upp2(), pett, N)  # for Julia 0.5

  push!(parentStack, startett.cliq )

  prepDwnPreOrderStack!(startett.bt, parentStack)
  processPreOrderStack!(startett.fg, startett.bt, parentStack, refdict, N=N, dbg=dbg, drawpdf=drawpdf )

  @info "dwnward leftovers, $(keys(refdict))"
  nothing
end

function prepPostOrderUpPassStacks!(bt::AbstractBayesTree,
                                    parentStack::Array{TreeClique,1},
                                    childStack::Array{TreeClique,1}  )
  # upward message passing preparation
  while ( length(parentStack) != 0 )
      #2.1 Pop a node from first stack and push it to second stack
      cliq = parentStack[end]
      deleteat!(parentStack, length(parentStack))
      push!(childStack, cliq )

      #2.2 Push left and right children of the popped node to first stack
      for child in out_neighbors(cliq, bt.bt) # nodedata.cliq, nodedata.bt.bt
          push!(parentStack, child )
      end
  end
  for child in childStack
    @show getLabel(child)
  end
  nothing
end

"""
    $SIGNATURES

Asynchronously perform up message passing, based on previoulsy prepared `chldstk::Vector{TreeClique}`.
"""
function asyncProcessPostStacks!(fgl::G,
                                 bt::AbstractBayesTree,
                                 chldstk::Vector{TreeClique},
                                 stkcnt::Int,
                                 refdict::Dict{Int,Future};
                                 N::Int=100,
                                 dbg::Bool=false,
                                 drawpdf::Bool=false  ) where G <: AbstractDFG
  #
  if stkcnt == 0
    @info "asyncProcessPostStacks! ERROR stkcnt=0"
    error("asyncProcessPostStacks! stkcnt=0")
  end
  cliq = chldstk[stkcnt]
  gomulti = true
  @info "Start Clique $(getLabel(cliq)) ============================="
  childMsgs = Array{NBPMessage,1}()
  ur = nothing
  for child in out_neighbors(cliq, bt.bt)
      @info "asyncProcessPostStacks -- $(stkcnt), cliq=$(getLabel(cliq)), start on child $(getLabel(child)) haskey=$(haskey(child.attributes, "remoteref"))"
        while !haskey(refdict, child.index)
          # info("Sleeping $(getLabel(cliq)) on lack of remoteref from $(getLabel(child))")
          # @show child.index, keys(refdict)
          sleep(0.25)
        end

      if gomulti
        ur = fetch(refdict[child.index])
      else
        ur = child.attributes["remoteref"]
      end
      updateFGBT!( fgl, bt, child.index, ur, dbg=dbg, fillcolor="pink" ) # deep copies happen in the update function
      drawpdf ? drawTree(bt) : nothing
      #delete!(child.attributes, "remoteref")

      push!(childMsgs, ur.upMsgs)

  end
  @info "====================== Clique $(getLabel(cliq)) ============================="
  emptr = BayesTree(nothing, 0, Dict{Int,TreeClique}(), Dict{String,Int}());
  pett = partialExploreTreeType(fgl, emptr, cliq, nothing, childMsgs) # bt   # parent cliq pointer is not needed here, fix Graphs.jl first

  if haskey(cliq.attributes, "remoteref")
      @info "asyncProcessPostStacks! -- WHY YOU ALREADY HAVE REMOTEREF?"
  end

  newprocid = upp2()
  if gomulti
    refdict[cliq.index] = remotecall(upGibbsCliqueDensity, newprocid, pett, N, dbg ) # swap order for Julia 0.5
  else
    # bad way to do non multi test
    cliq.attributes["remoteref"] = upGibbsCliqueDensity(pett, N, dbg)
  end

  # delete as late as possible, but could happen sooner
  for child in out_neighbors(cliq, bt.bt)
      # delete!(child.attributes, "remoteref")
      delete!(refdict, child.index)
  end

  @info "End Clique $(getLabel(cliq)) ============================="
  nothing
end

"""
    $SIGNATURES

Multiprocess upward belief propagation message passing function, using async tasks.

Notes
- asyncs used to wrap remotecall for multicore.
- separate multithreaded calls can occur on each separate process.
"""
function processPostOrderStacks!(fg::G,
                                 bt::AbstractBayesTree,
                                 childStack::Array{TreeClique,1};
                                 N::Int=100,
                                 dbg::Bool=false,
                                 drawpdf::Bool=false  ) where G <: AbstractDFG
  #
  refdict = Dict{Int,Future}()

  stkcnt = length(childStack)
  @sync begin
    sendcnt = stkcnt:-1:1 # separate stable memory
    for i in 1:stkcnt
        @async asyncProcessPostStacks!(fg, bt, childStack, sendcnt[i], refdict, N=N, dbg=dbg, drawpdf=drawpdf ) # deepcopy(stkcnt)
    end
  end
  @info "processPostOrderStacks! -- THIS ONLY HAPPENS AFTER SYNC"
  # we still need to fetch the root node computational output
  if true
    ur = fetch(refdict[childStack[1].index])
  else
    ur = childStack[1].attributes["remoteref"]
  end
  # delete!(childStack[1].attributes, "remoteref") # childStack[1].cliq
  delete!(refdict, childStack[1].index)

  @info "upward leftovers, $(keys(refdict))"

  updateFGBT!(fg, bt, childStack[1].index, ur, dbg=dbg, fillcolor="pink" ) # nodedata
  drawpdf ? drawTree(bt) : nothing
  nothing
end


"""
    $SIGNATURES

Return clique pointers for the given order in which they will be solved (sequentially).
"""
function getCliqOrderUpSolve(treel::AbstractBayesTree, startcliq=treel.cliques[1])
  # http://www.geeksforgeeks.org/iterative-postorder-traversal/
  # this is where we launch the downward iteration process from
  parentStack = Vector{TreeClique}()
  childStack = Vector{TreeClique}()
  #Loop while first stack is not empty
  push!(parentStack, startcliq)
  # Starting at the root means we have a top down view of the tree
  prepPostOrderUpPassStacks!(treel, parentStack, childStack)
  return childStack
end

"""
    $SIGNATURES

Return clique pointers for the given order in which they will be solved (sequentially).
"""
getTreeCliqSolveOrderUp(treel::AbstractBayesTree, startcliq=treel.cliques[1]) = getCliqOrderUpSolve(treel, startcliq)

"""
    $SIGNATURES

Perform upward message passing (multi-process) algorithm for sum-product solution from leaves to root of the tree.

Notes:
* inspired by http://www.geeksforgeeks.org/iterative-postorder-traversal/
* this is where downward iteration process is launched from.
"""
function upMsgPassingIterative!(startett::ExploreTreeType{T};
                                N::Int=100,
                                dbg::Bool=false,
                                drawpdf::Bool=false  ) where {T}
  #
  childStack = getCliqOrderUpSolve(startett.bt, startett.cliq)
  # Starting at the root means we have a top down view of the tree
  processPostOrderStacks!(startett.fg, startett.bt, childStack, N=N, dbg=dbg, drawpdf=drawpdf)
  nothing
end
# for (ids, cliq) in treel.cliques
#   getData(cliq).initialized = :initialized
# end




"""
    $SIGNATURES

List vertices two neighbors deep.
"""
function ls2(fgl::FactorGraph, vsym::Symbol)
  @warn "Deprecated, please use DFG.ls2"
  xxf = ls(fgl, vsym)
  xlxl = Symbol[]
  for xf in xxf
    xx = lsf(fgl,xf)
    xlxl = union(xlxl, xx)
  end
  xlxl = setdiff(xlxl, [vsym])
  return xlxl
end
