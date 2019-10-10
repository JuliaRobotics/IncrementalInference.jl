
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














# explore all shortest paths combinations in verts, add neighbors and reference subgraph
# Using unique index into graph data structure
function subgraphFromVerts(fgl::FactorGraph,
    verts::Array{Int,1};
    neighbors::Int=0  )

  vertdict = Dict{Int,Graphs.ExVertex}()
  for vert in verts
    vertdict[vert] = fgl.g.vertices[vert]
  end

  return subgraphFromVerts(fgl,vertdict,neighbors=neighbors)
end

"""
    $(SIGNATURES)

Explore all shortest paths combinations in verts, add neighbors and reference
subgraph using unique index into graph data structure.
"""
function subGraphFromVerts(fgl::FactorGraph,
                           verts::T;
                           neighbors::Int=0  ) where {T <: Union{Vector{String},Vector{Symbol}}}
  #
  vertdict = Dict{Int,Graphs.ExVertex}()
  for vert in verts
    id = -1
    vsym = Symbol(vert)
    if haskey(fgl.IDs, vsym)
      id = fgl.IDs[vsym]
    else
      error("FactorGraph01 only supports variable node subgraph search at this point")
    end
    vertdict[id] = getVert(fgl, vsym) # fgl.g.vertices[id]
  end

  return subgraphFromVerts(fgl,vertdict,neighbors=neighbors)
end


function subgraphFromVerts(fgl::FactorGraph,
                           verts::T;
                           neighbors::Int=0  ) where {T <: Union{Vector{String},Vector{Symbol}, Dict{Int,Graphs.ExVertex}}}
  #
  @warn "`subgraphFromVerts` deprecated, use `subGraphFromVerts` instead."
  subGraphFromVerts(fgl, verts, neighbors=neighbors)
end






"""
    $SIGNATURES

Return array of all variable vertices in a clique.
"""
function getCliqVars(subfg::FactorGraph, cliq::Graphs.ExVertex)
  verts = Graphs.ExVertex[]
  for vid in getCliqVars(subfg, cliq)
    push!(verts, getVert(subfg, vid, api=localapi))
  end
  return verts
end




### TODO: TO BE REFACTORED FOR DFG

# some plotting functions on the factor graph
function stackVertXY(fg::FactorGraph, lbl::String)
    v = dlapi.getvertex(fg,lbl)
    vals = getVal(v)
    X=vec(vals[1,:])
    Y=vec(vals[2,:])
    return X,Y
end

### TODO: TO BE REFACTORED FOR DFG

"""
    $(SIGNATURES)

Save mostly complete Factor Graph type by converting complicated FunctionNodeData
types to 'Packed' types using user supplied converters. Ground truth can also be
saved and recovered by the associated loadjld(file="tempfg.jld2") method.

Notes:
- Must use `.jld2` since Julia 1.0 (previous version was deprecated).
"""
function savejld(fgl::G;
                 file::AbstractString="tempfg.jld2"  ) where G <: AbstractDFG
  #
  @error "savejld has been deprecated, use saveDFG instead"
  # fgs = encodefg(fgl)
  @save file fgl
  return file
end

"""
    $(SIGNATURES)

Opposite of savejld(fg, gt=gt, file="tempfg.jl") to load data from file. This function
uses the unpacking converters for converting all PackedInferenceType to FunctorInferenceType.
"""
function loadjld(;file::AbstractString="tempfg.jld2")
  @error "loadjld has been deprecated, use loadDFG instead"
  fgd = @load file fgl
  return fgd
end


"""
    $(SIGNATURES)

Return the last up message stored in `cliq` of Bayes (Junction) tree.
"""
function upMsg(cliq::Graphs.ExVertex)
  @warn "deprecated upMsg, use getUpMsg instead"
  getData(cliq).upMsg
end
function upMsg(btl::BayesTree, sym::Symbol)
  @warn "deprecated upMsg, use getUpMsg instead"
  upMsg(whichCliq(btl, sym))
end

"""
    $(SIGNATURES)

Return the last down message stored in `cliq` of Bayes (Junction) tree.
"""
function dwnMsg(cliq::Graphs.ExVertex)
  @warn "deprecated dwnMsg, use getDwnMsgs instead"
  getData(cliq).dwnMsg
end
function dwnMsg(btl::BayesTree, sym::Symbol)
  @warn "deprecated dwnMsg, use getDwnMsgs instead"
  dwnMsg(whichCliq(btl, sym))
end


function getCliquePotentials!(fg::FactorGraph, bt::BayesTree, chkcliq::Int)
    @error "getCliquePotentials! deprecated, use setCliqPotentials! with DFG objects instead of FactorGraph"
    setCliqPotentials!(fg, bt.cliques[chkcliq])
end

"""
    $(SIGNATURES)

Pass NBPMessages back down the tree -- pre order tree traversal.
"""
function downMsgPassingRecursive(inp::ExploreTreeType{T}; N::Int=100, dbg::Bool=false, drawpdf::Bool=false) where {T}
  @info "====================== Clique $(inp.cliq.attributes["label"]) ============================="

  mcmciter = inp.prnt != nothing ? 3 : 0; # skip mcmc in root on dwn pass
  rDDT = downGibbsCliqueDensity(inp.fg, inp.cliq, inp.sendmsgs, N, mcmciter, dbg) #dwnMsg
  updateFGBT!(inp.fg, inp.bt, inp.cliq.index, rDDT, dbg=dbg, fillcolor="lightblue")
  setCliqStatus!(inp.cliq, :downsolved)
  drawpdf ? drawTree(inp.bt) : nothing

  # rr = Array{Future,1}()
  # pcs = procs()

  ddt=nothing
  for child in out_neighbors(inp.cliq, inp.bt.bt)
    ett = ExploreTreeType(inp.fg, inp.bt, child, inp.cliq, [rDDT.dwnMsg])#inp.fg
    ddt = downMsgPassingRecursive( ett , N=N, dbg=dbg )
    drawpdf ? drawTree(inp.bt) : nothing
  end

  # return modifications to factorgraph to calling process
  return ddt
end

# post order tree traversal and build potential functions
function upMsgPassingRecursive(inp::ExploreTreeType{T}; N::Int=100, dbg::Bool=false, drawpdf::Bool=false) where {T}
    @info "Start Clique $(inp.cliq.attributes["label"]) ============================="
    childMsgs = Array{NBPMessage,1}()

    outnei = out_neighbors(inp.cliq, inp.bt.bt)
    len = length(outnei)
    for child in outnei
        ett = ExploreTreeType(inp.fg, inp.bt, child, inp.cliq, NBPMessage[])
        @info "upMsgRec -- calling new recursive on $(ett.cliq.attributes["label"])"
        newmsgs = upMsgPassingRecursive(  ett, N=N, dbg=dbg ) # newmsgs
        @info "upMsgRec -- finished with $(ett.cliq.attributes["label"]), w $(keys(newmsgs.p)))"
        push!(  childMsgs, newmsgs )
    end

    @info "====================== Clique $(inp.cliq.attributes["label"]) ============================="
    ett = ExploreTreeType(inp.fg, inp.bt, inp.cliq, nothing, childMsgs)

    urt = upGibbsCliqueDensity(ett, N, dbg) # upmsgdict
    updateFGBT!(inp.fg, inp.bt, inp.cliq.index, urt, dbg=dbg, fillcolor="lightblue")
    drawpdf ? drawTree(inp.bt) : nothing
    @info "End Clique $(inp.cliq.attributes["label"]) ============================="
    urt.upMsgs
end


function downGibbsCliqueDensity(inp::ExploreTreeType{T},
                                N::Int=100,
                                dbg::Bool=false,
                                logger=ConsoleLogger()  ) where {T}
  #
  with_logger(logger) do
    @info "=================== Iter Clique $(inp.cliq.attributes["label"]) ==========================="
  end
  mcmciter = inp.prnt != nothing ? 3 : 0
  return downGibbsCliqueDensity(inp.fg, inp.cliq, inp.sendmsgs, N, mcmciter, dbg)
end

function prepDwnPreOrderStack!(bt::BayesTree,
                               parentStack::Array{Graphs.ExVertex,1})
  # dwn message passing function
  nodedata = nothing
  tempStack = Array{Graphs.ExVertex,1}()
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

function findVertsAssocCliq(fgl::FactorGraph, cliq::Graphs.ExVertex)

  cliqdata = getData(cliq)
  IDS = [cliqdata.frontalIDs; cliqdata.conditIDs] #inp.cliq.attributes["frontalIDs"]


  @error "findVertsAssocCliq -- not completed yet"
  nothing
end

function partialExploreTreeType(pfg::G, pbt::BayesTree, cliqCursor::Graphs.ExVertex, prnt, pmsgs::Array{NBPMessage,1}) where G <: AbstractDFG
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
                             bt::BayesTree,
                             parentStack::Array{Graphs.ExVertex,1},
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

  emptr = BayesTree(nothing, 0, Dict{Int,Graphs.ExVertex}(), Dict{String,Int}());

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
                               bt::BayesTree,
                               parentStack::Array{Graphs.ExVertex,1},
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
  parentStack = Array{Graphs.ExVertex,1}()
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

function prepPostOrderUpPassStacks!(bt::BayesTree,
                                    parentStack::Array{Graphs.ExVertex,1},
                                    childStack::Array{Graphs.ExVertex,1}  )
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
    @show child.attributes["label"]
  end
  nothing
end

"""
    $SIGNATURES

Asynchronously perform up message passing, based on previoulsy prepared `chldstk::Vector{ExVertex}`.
"""
function asyncProcessPostStacks!(fgl::G,
                                 bt::BayesTree,
                                 chldstk::Vector{Graphs.ExVertex},
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
  @info "Start Clique $(cliq.attributes["label"]) ============================="
  childMsgs = Array{NBPMessage,1}()
  ur = nothing
  for child in out_neighbors(cliq, bt.bt)
      @info "asyncProcessPostStacks -- $(stkcnt), cliq=$(cliq.attributes["label"]), start on child $(child.attributes["label"]) haskey=$(haskey(child.attributes, "remoteref"))"
        while !haskey(refdict, child.index)
          # info("Sleeping $(cliq.attributes["label"]) on lack of remoteref from $(child.attributes["label"])")
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
  @info "====================== Clique $(cliq.attributes["label"]) ============================="
  emptr = BayesTree(nothing, 0, Dict{Int,Graphs.ExVertex}(), Dict{String,Int}());
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

  @info "End Clique $(cliq.attributes["label"]) ============================="
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
                                 bt::BayesTree,
                                 childStack::Array{Graphs.ExVertex,1};
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
function getCliqOrderUpSolve(treel::BayesTree, startcliq=treel.cliques[1])
  # http://www.geeksforgeeks.org/iterative-postorder-traversal/
  # this is where we launch the downward iteration process from
  parentStack = Vector{Graphs.ExVertex}()
  childStack = Vector{Graphs.ExVertex}()
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
getTreeCliqSolveOrderUp(treel::BayesTree, startcliq=treel.cliques[1]) = getCliqOrderUpSolve(treel, startcliq)

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
$(SIGNATURES)

Add a node (variable) to a graph. Use this over the other dispatches.

DEPRECATED: use addVarialbe! instead.
"""
function addNode!(fg::FactorGraph,
                  lbl::Symbol,
                  softtype::T;
                  N::Int=100,
                  autoinit::Bool=true,  # does init need to be separate from ready? TODO
                  ready::Int=1,
                  dontmargin::Bool=false,
                  labels::Vector{<:AbstractString}=String[],
                  uid::Int=-1,
                  smalldata=""  ) where {T <:InferenceVariable}
  #
  @warn "IIF.addNode!(..) is being deprecated, use IIF.addVariable!(..) instead."
  return addVariable!( fg,
                       lbl,
                       softtype,
                       N=N,
                       autoinit=autoinit,  # does init need to be separate from ready? TODO
                       ready=ready,
                       dontmargin=dontmargin,
                       labels=labels,
                       uid=uid,
                       smalldata=smalldata )
end
function addNode!(fg::FactorGraph,
                  lbl::Symbol,
                  softtype::Type{<:InferenceVariable};
                  N::Int=100,
                  autoinit::Bool=true,
                  ready::Int=1,
                  dontmargin::Bool=false,
                  labels::Vector{<:AbstractString}=String[],
                  uid::Int=-1,
                  # dims::Int=-1,
                  smalldata=""  )
  #
  @warn "IIF.addNode!(..) is being deprecated, use IIF.addVariable!(..) instead."
  return addVariable!(fg,
                      lbl,
                      softtype,
                      N=N,
                      autoinit=autoinit,
                      ready=ready,
                      dontmargin=dontmargin,
                      labels=labels,
                      uid=uid,
                      smalldata=smalldata  )
end



"""
    $(SIGNATURES)

Unpack PackedFunctionNodeData formats back to regular FunctonNodeData.
"""
function decodefg(fgs::FactorGraph)
  fgu = deepcopy(fgs)
  fgu.cg = nothing # will be deprecated or replaced
  fgu.registeredModuleFunctions = nothing # TODO: obsolete
  fgu.g = Graphs.incdict(Graphs.ExVertex,is_directed=false)
  @showprogress 1 "Decoding variables..." for (vsym,vid) in fgs.IDs
    cpvert = deepcopy(getVert(fgs, vid, api=api))
    api.addvertex!(fgu, cpvert) #, labels=vnlbls)  # currently losing labels
  end

  @showprogress 1 "Decoding factors..." for (fsym,fid) in fgu.fIDs
    fdata = solverData(fgs, fid)
    data = decodePackedType(fdata, "")

    # data = FunctionNodeData{ftyp}(Int[], false, false, Int[], m, gwpf)
    newvert = ExVertex(fid,string(fsym))
    for (key,val) in getVert(fgs,fid,api=api).attributes
      newvert.attributes[key] = val
    end
    setData!(newvert, data)
    api.addvertex!(fgu, newvert)
  end
  fgu.g.inclist = typeof(fgs.g.inclist)()

  # iterated over all edges
  @showprogress 1 "Decoding edges..." for (eid, edges) in fgs.g.inclist
    fgu.g.inclist[eid] = Vector{typeof(edges[1])}()
    for ed in edges
      newed = Graphs.Edge(ed.index,
          fgu.g.vertices[ed.source.index],
          fgu.g.vertices[ed.target.index]  )
      push!(fgu.g.inclist[eid], newed)
    end
  end

  # rebuild factormetadata
  @showprogress 1 "Rebuilding factor metadata..." for (fsym,fid) in fgu.fIDs
    varuserdata = []
    fcnode = getVert(fgu, fsym, nt=:fnc)
    # ccw = solverData(fcnode)
    ccw_jld = deepcopy(solverData(fcnode))
    allnei = Graphs.ExVertex[]
    for nei in out_neighbors(fcnode, fgu.g)
        push!(allnei, nei)
        data = IncrementalInference.solverData(nei)
        push!(varuserdata, data.softtype)
    end
    setDefaultFactorNode!(fgu, fcnode, allnei, ccw_jld.fnc.usrfnc!, threadmodel=ccw_jld.fnc.threadmodel, multihypo=veeCategorical(ccw_jld.fnc.hypotheses))
    ccw_new = IncrementalInference.solverData(fcnode)
    for i in 1:Threads.nthreads()
      ccw_new.fnc.cpt[i].factormetadata.variableuserdata = deepcopy(varuserdata)
    end
    ## Rebuild solverData(fcnode).fncargvID, however, the list is order sensitive
    # out_neighbors does not gaurantee ordering -- i.e. why is it not being saved
    for field in fieldnames(typeof(ccw_jld))
      if field != :fnc
        setfield!(ccw_new, field, getfield(ccw_jld, field))
      end
    end
  end
  return fgu
end



# """
#     $(SIGNATURES)
#
# Return all elements `ls(fg)` as tuples, or nodes connected to the a specific element, eg. `ls(fg, :x1)
# """
# function ls(fgl::FactorGraph, lbl::Symbol; ring::Int=1)
#   @warn "Deprecated, please use DFG.ls"
#   # TODO ring functionality must still be implemented
#   lsa = Symbol[]
#   # v = nothing
#   if haskey(fgl.IDs, lbl)
#     id = fgl.IDs[lbl]
#   else
#     return lsa
#   end
#   # this is unnecessary
#   v = getVariable(fgl,id)
#   for outn in api.outneighbors(fgl, v)
#     # if outn.attributes["ready"] = 1 && outn.attributes["backendset"]=1
#       push!(lsa, Symbol(outn.label))
#     # end
#   end
#   return lsa
# end
# ls(fgl::FactorGraph, lbl::T) where {T <: AbstractString} = ls(fgl, Symbol(lbl))
#
# """
#     $(SIGNATURES)
#
# Experimental union of elements version of ls(::FactorGraph, ::Symbol).  Not mean't to replace broadcasting `ls.(fg, [:x1;:x2])`
# """
# function ls(fgl::FactorGraph,
#             lbls::Vector{Symbol};
#             ring::Int=1)
#   @warn "Deprecated, please use DFG.ls"
#   union(ls.(fgl, lbls, ring=ring)[:]...)
# end
#
# """
#     $(SIGNATURES)
#
# List the nodes in a factor graph.
#
# # Examples
# ```julia-repl
# ls(fg)
# ```
# """
# function ls(fgl::FactorGraph; key1='x', key2='l')
#   @warn "Deprecated, please use DFG.ls"
#   k = collect(keys(fgl.IDs))
#   x, l = String[], String[]
#   xval, lval = Int[], Int[]
#   xstr, lstr = String[], String[]
#   xvalnested, lvalnested = String[], String[]
#   xstrnested, lstrnested = String[], String[]
#   canparse1, canparse2 = true,true
#   nestedparse1, nestedparse2 = true, true
#   idx = 0
#   for id in k
#     idx += 1
#     idstr = string(id)
#     # val = parse(Int,kstr[2:end]) # TODO: handle non-int labels
#     node_idx = idstr[2:end]
#     canparse = allnums(node_idx)
#     nested = isnestednum(node_idx)
#     if idstr[1] == key1
#       keystr = string(key1,node_idx)
#       if canparse
#         push!(xstr, keystr)
#         push!(xval, parse(Int, node_idx))
#       elseif nested
#         push!(xvalnested, node_idx)
#         push!(xstrnested, string(node_idx))
#       else
#         push!(x,keystr)
#       end
#     elseif idstr[1] == key2
#       keystr = string(key2,node_idx)
#       if canparse
#         push!(lstr, keystr)
#         push!(lval, parse(Int, node_idx))
#       elseif nested
#         push!(lstrnested, keystr)
#         push!(lvalnested, string(node_idx))
#       else
#         push!(l,string(key2,node_idx))
#       end
#     end
#   end
#   x1 = xstr[sortperm(xval)]
#   x2 = xstrnested[sortnestedperm(xvalnested)]
#   x = [x1; x2; sort(x)]
#
#   l1 = lstr[sortperm(lval)]
#   l2 = lstrnested[sortnestedperm(lvalnested)]
#   l = [l1; l2; sort(l)]
#
#   xx = Symbol.(x)
#   ll = Symbol.(l)
#   return xx, ll #return poses, landmarks
# end
#
# lsf(fgl::FactorGraph) = collect(keys(fgl.fIDs))
#
# """
#     $(SIGNATURES)
#
# List factors in a factor graph.
#
# # Examples
# ```julia-repl
# lsf(fg, :x1)
# ```
# """
# function lsf(fgl::FactorGraph, lbl::Symbol)
#   @warn "Deprecated, please use DFG.lsf"
#   lsa = Symbol[]
#   if haskey(fgl.fIDs, lbl)
#     id = fgl.fIDs[lbl]
#   else
#     return lsa
#   end
#   v = getVariable(fgl, id) # fgl.g.vertices[id] #fgl.f[id]
#   for outn in api.outneighbors(fgl, v) # out_neighbors(v, fgl.g)
#     push!(lsa, Symbol(outn.label))
#   end
#   return lsa
# end
#
#
# """
#     $(SIGNATURES)
#
# List factors in a factor graph.
#
# # Examples
# ```julia-repl
# lsf(fg)
# ```
# """
# lsf(fgl::FactorGraph, lbl::T) where {T <: AbstractString} = lsf(fgl,Symbol(lbl))
#
# function lsf(fgl::FactorGraph,
#              mt::Type{T}) where {T <: FunctorInferenceType}
#   @warn "Deprecated, please use DFG.lsf"
#   syms = Symbol[]
#   for (fsym,fid) in fgl.fIDs
#     if typeof(getfnctype(fgl, fid))==T
#       push!(syms, fsym)
#     end
#   end
#   return syms
# end
#
# function lsf(fgl::FactorGraph)
#   @warn "Deprecated, please use DFG.lsf"
#   collect(keys(fgl.fIDs))
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

function initializeNode!(fgl::G,
                         sym::Symbol;
                         N::Int=100  ) where G <: AbstractDFG
  #
  @warn "initializeNode! has been deprecated in favor of initVariable!"
  initVariable!(fgl,sym,N=N )
end



# """
#     $(SIGNATURES)
#
# Test if all elements of the string is a number:  Ex, "123" is true, "1_2" is false.
# """
# allnums(str::S) where {S <: AbstractString} = occursin(Regex(string(["[0-9]" for j in 1:length(str)]...)), str)
# # occursin(r"_+|,+|-+", node_idx)
#
# isnestednum(str::S; delim='_') where {S <: AbstractString} = occursin(Regex("[0-9]+$(delim)[0-9]+"), str)
#
# function sortnestedperm(strs::Vector{<:AbstractString}; delim='_')
#   str12 = split.(strs, delim)
#   sp1 = sortperm(parse.(Int,getindex.(str12,2)))
#   sp2 = sortperm(parse.(Int,getindex.(str12,1)[sp1]))
#   return sp1[sp2]
# end

# moved to DFG
# """
#     $SIGNATURES
#
# Return `::Bool` on whether `fg::FactorGraph` has orphaned nodes or graph fragments.
# """
# hasOrphans(fg::FactorGraph) = sum(length.(ls.(fg, [ls(fg)[1];ls(fg)[2]])) .== 0) > 0
