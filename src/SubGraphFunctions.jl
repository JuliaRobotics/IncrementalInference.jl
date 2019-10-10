
# TODO -- convert to use add_vertex! instead, since edges type must be made also
function addVerticesSubgraph(fgl::FactorGraph,
    fgseg::FactorGraph,
    vertdict::Dict{Int,Graphs.ExVertex})

    for vert in vertdict
      fgseg.g.vertices[vert[1]] = vert[2]
      if haskey(fgl.v,vert[1])
        fgseg.g.vertices[vert[1]] = vert[2]
        fgseg.IDs[Symbol(vert[2].label)] = vert[1]

        # add edges going in opposite direction
        elr = Graphs.out_edges(vert[2], fgl.g)
        len = length(elr)
        keeprm = trues(len)
        j = 0
        for i in 1:len
          if !haskey(vertdict, elr[i].target.index) # a function node in set, so keep ref
            keeprm[i] = false
            j+=1
          end
        end
        if j < len
          elridx = elr[1].source.index
          fgseg.g.inclist[elridx] = elr[keeprm]
        end
      elseif haskey(fgl.f, vert[1])
        fgseg.f[vert[1]] = vert[2] # adding element to subgraph
        fgseg.fIDs[Symbol(vert[2].label)] = vert[1]
        # get edges associated with function nodes and push edges onto incidence list
        el = Graphs.out_edges(vert[2], fgl.g)
        elidx = el[1].source.index
        fgseg.g.inclist[elidx] = el # okay because treating function nodes only
        fgseg.g.nedges += length(el)
      else
        error("Unknown type factor graph vertex type, something is wrong")
      end
    end
    nothing
end

"""
    $SIGNATURES

Construct a new factor graph object as a subgraph of `fgl::FactorGraph` based on the
variable labels `syms::Vector{Symbols}`.

Notes
- Slighly messy internals, but gets the job done -- some room for performance improvement.

Related

getVariableIds
"""
function buildSubgraphFromLabels(dfg::G,
                                 syms::Vector{Symbol},
                                 destType::Type{<:AbstractDFG}=InMemDFGType ) where G <: AbstractDFG

  # data structure for cliq sub graph
  if G <: InMemoryDFGTypes
    #Same type
    cliqSubFg = initfg(G)
  else
    #Default
    cliqSubFg = initfg(destType)
  end

  # add a little too many variables (since we need the factors)
  for sym in syms
    DFG.getSubgraphAroundNode(dfg, DFG.getVariable(dfg, sym), 2, false, cliqSubFg)
  end

  # remove excessive variables that were copied by neighbors distance 2
  currVars = DFG.getVariableIds(cliqSubFg)
  toDelVars = setdiff(currVars, syms)
  for dv in toDelVars
    # delete any neighboring factors first
    for fc in DFG.lsf(cliqSubFg, dv)
      DFG.deleteFactor!(cliqSubFg, fc)
    end

    # and the variable itself
    DFG.deleteVariable!(cliqSubFg, dv)
  end

  return cliqSubFg
end
# function buildSubgraphFromLabels(dfg::G, syms::Vector{Symbol}) where G <: AbstractDFG
#   fgseg = initfg() #sessionname=dfg.sessionname, robotname=dfg.robotname)
#
#   for sym in syms
#     vert = DFG.getVariable(dfg, sym) #, api=localapi)
#     st = getSofttype(vert)
#     addVariable!(fgseg, sym, st) #, uid=vert.index)
#     if isInitialized(dfg,sym)
#       manualinit!(fgseg, sym, getKDE(vert))
#     end
#   end
#
#   for sym in syms
#     for fct in DFG.ls(dfg, :x1)
#       if !hasFactor(fgseg, fct)
#         # check all variables are in desired variable set
#         possibleVars = DFG.lsf(dfg, fct)
#         ivars = intersect(possibleVars, syms)
#         @show length(ivars), length(possibleVars)
#         if length(ivars) == length(possibleVars)
#           # fvert = getVert(dfg, fct, api=localapi, nt=:fct)
#           ufc = DFG.getFactor(dfg, fct) # fvert
#
#           addFactor!(fgseg, possibleVars, getData(ufc).fnc.usrfnc!, autoinit=false) #, uid=fvert.index)
#         end
#       end
#     end
#   end
#
#   return fgseg
# end



# NOTICE, nodeIDs and factorIDs are not brough over by this method yet
# must sort out for incremental updates
function genSubgraph(fgl::FactorGraph,
                     vertdict::Dict{Int,Graphs.ExVertex})
  # edgedict::Dict{Int,Graphs.Edge{Graphs.ExVertex}},
  #

  fgseg = initfg(sessionname=fgl.sessionname, robotname=fgl.robotname)

  error("genSubgraph() is obsolete")

  # fgseg = FactorGraph() # new handle for just a segment of the graph
  # fgseg.g = Graphs.inclist(Graphs.ExVertex, is_directed=false)
  #
  # fgseg.v = Dict{Int,Graphs.ExVertex}()
  # fgseg.f = Dict{Int,Graphs.ExVertex}()
  # fgseg.IDs = Dict{AbstractString,Int}()
  # fgseg.fIDs = Dict{AbstractString,Int}()
  #
  # # TODO -- convert to use empty constructor since Graphs.incdict now works
  # fgseg.g.vertices = Array{Graphs.ExVertex,1}(length(fgl.g.vertices))
  # fgseg.g.inclist = Array{Array{Graphs.Edge{Graphs.ExVertex},1},1}(length(fgl.g.inclist))
  #
  # addVerticesSubgraph(fgl, fgseg, vertdict)
  #
  # fgseg.id = fgl.id
  # fgseg.bnid = fgl.bnid
  # fgseg.dimID = fgl.dimID

  return fgseg
end

function getShortestPathNeighbors(fgl::FactorGraph;
    from::Graphs.ExVertex=nothing,
    to::Graphs.ExVertex=nothing,
    neighbors::Int=0 )

  edgelist = shortest_path(fgl.g, ones(num_edges(fgl.g)), from, to)
  vertdict = Dict{Int,Graphs.ExVertex}()
  edgedict = edgelist2edgedict(edgelist)
  expandVertexList!(fgl, edgedict, vertdict) # grow verts
  for i in 1:neighbors
    expandEdgeListNeigh!(fgl, vertdict, edgedict) # grow edges
    expandVertexList!(fgl, edgedict, vertdict) # grow verts
  end
  return vertdict
end

function subgraphShortestPath(fgl::FactorGraph;
                              from::Graphs.ExVertex=nothing,
                              to::Graphs.ExVertex=nothing,
                              neighbors::Int=0  )
  #
  vertdict = getShortestPathNeighbors(fgl, from=from, to=to, neighbors=neighbors)
  return genSubgraph(fgl, vertdict)
end

# explore all shortest paths combinations in verts, add neighbors and reference subgraph
function subGraphFromVerts(fgl::FactorGraph,
                           verts::Dict{Int,Graphs.ExVertex};
                           neighbors::Int=0  )
  #
  allverts = Dict{Int,Graphs.ExVertex}()
  allkeys = collect(keys(verts))
  len = length(allkeys)
  # union all shortest path combinations in a vertdict
  for i in 1:len, j in (i+1):len
    from = verts[allkeys[i]]
    to = verts[allkeys[j]]
    vertdict = getShortestPathNeighbors(fgl, from=from, to=to, neighbors=neighbors)
    for vert in vertdict
      if !haskey(allverts, vert[1])
        allverts[vert[1]] = vert[2]
      end
    end
  end

  return genSubgraph(fgl, allverts)
end




"""
    $SIGNATURES

Transfer contents of `src::FactorGraph` variables `syms::Vector{Symbol}` to `dest::FactorGraph`.

Notes
- Reads, `dest` := `src`, for all `syms`
"""
function transferUpdateSubGraph!(dest::G1,
                                 src::G2,
                                 syms::Vector{Symbol}=union(ls(src)...),
                                 logger=ConsoleLogger()  ) where {G1 <: AbstractDFG, G2 <: AbstractDFG}
  #
  with_logger(logger) do
    @info "transferUpdateSubGraph! -- syms=$syms"

    # TODO add with DFG v0.4
    # DFG.updateGraphSolverData!(src, dest, syms)
    for sym in syms
      vari = DFG.getVariable(src, sym)
      rc = size(solverData(vari).val)
      pp = getKDE(vari)
      rc2 = size(getPoints(pp))
      @info "sym=$sym, mem size of val=$rc and $(rc2)"
      updateFullVertData!(dest, vari, updateMAPest=true)
    end
  end
  nothing
end
