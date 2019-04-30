
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
"""
function buildSubgraphFromLabels(fgl::FactorGraph, syms::Vector{Symbol})
  fgseg = initfg(sessionname=fgl.sessionname, robotname=fgl.robotname)

  for sym in syms
    vert = getVert(fgl, sym, api=localapi)
    st = getSofttype(vert)
    addVariable!(fgseg, sym, st, uid=vert.index)
    manualinit!(fgseg, sym, getKDE(vert))
  end

  for sym in syms
    for fct in ls(fgl, sym)
      if !hasFactor(fgseg, fct)
        # check all variables are in desired variable set
        possibleVars = lsf(fgl, fct)
        ivars = intersect(possibleVars, syms)
        if length(ivars) == length(possibleVars)
          fvert = getVert(fgl, fct, api=localapi, nt=:fct)
          ufc = getFactor(fvert)
          addFactor!(fgseg, possibleVars, ufc, autoinit=false, uid=fvert.index)
        end
      end
    end
  end

  return fgseg
end



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
