
reshapeVec2Mat(vec::Vector, rows::Int) = reshape(vec, rows, round(Int,length(vec)/rows))
# function reshapeVec2Mat(vec::Vector, rows::Int)
#   M = reshape(vec, rows, round(Int,length(vec)/rows))
#   return ndims(M) < 2 ? (M')' : M
# end

# get vertex from factor graph according to label symbol "x1"
function getVert(fgl::FactorGraph, lbl::ASCIIString)
  # if haskey(fgl.IDs, lbl)
  #   return fgl.g.vertices[fgl.IDs[lbl]]
  # end
  dlapi.getvertex(fgl, lbl)
end
function getVert(fgl::FactorGraph, id::Int64)
  dlapi.getvertex(fgl, id)
end

getData(v::Graphs.ExVertex) = v.attributes["data"]

function getVal(v::Graphs.ExVertex)
  return getData(v).val
end

# Convenience function to get values for given variable label
function getVal(fgl::FactorGraph, lbl::ASCIIString)
  getVal(dlapi.getvertex(fgl, lbl))
end
function getVal(fgl::FactorGraph, exvertid::Int64)
  getVal(dlapi.getvertex(fgl, exvertid))
end


# setVal! assumes you will update values to database separate, this used for local graph mods only
function setVal!(v::Graphs.ExVertex, val::Array{Float64,2})
  v.attributes["data"].val = val
  # v.attributes["val"] = val
  nothing
end
function getBWVal(v::Graphs.ExVertex)
  return getData(v).bw
end
function setBW!(v::Graphs.ExVertex, bw::Array{Float64,2})
  v.attributes["data"].bw = bw
  # v.attributes["bw"] = bw
  nothing
end
function setVal!(v::Graphs.ExVertex, val::Array{Float64,2}, bw::Array{Float64,2})
  setVal!(v,val)
  setBW!(v,bw)
  nothing
end
function setVal!(v::Graphs.ExVertex, val::Array{Float64,2}, bw::Array{Float64,1})
  setVal!(v,val,(bw')')
  nothing
end
function setValKDE!(v::Graphs.ExVertex, val::Array{Float64,2})
  p = kde!(val)
  setVal!(v,val,getBW(p)[:,1]) # TODO -- this can be little faster
  nothing
end

# TODO -- there should be a better way, without retrieving full vertex
getOutNeighbors(fgl::FactorGraph, v::ExVertex) = dlapi.outneighbors(fgl,v)
getOutNeighbors(fgl::FactorGraph, vertid::Int64) = dlapi.outneighbors(fgl, dlapi.getvertex(fgl,vertid) )

function updateFullVert!(fgl::FactorGraph, exvert::ExVertex)
  dlapi.updatevertex!(fgl, exvert)
end


function setDefaultNodeData!(v::Graphs.ExVertex, initval::Array{Float64,2},
                              stdev::Array{Float64,2}, dodims::Int64, N::Int64)
  pN = Union{}
  if size(initval,2) < N
    p = kde!(initval,diag(stdev));
    pN = resample(p,N)
  else
    pN = kde!(initval, "lcv")
  end
  dims = size(initval,1) # rows indicate dimensions
  sp = round(Int64,linspace(dodims,dodims+dims-1,dims))
  data = VariableNodeData(initval, stdev, getPoints(pN),
                          (getBW(pN)[:,1]')', Int64[], sp,
                          dims, false, 0, Int64[])

  v.attributes["data"] = data

  nothing
end

function addNewVarVertInGraph!(fgl::FactorGraph, vert::Graphs.ExVertex, id::Int64, lbl::AbstractString, ready::Int)
  vert.attributes = Graphs.AttributeDict() #fg.v[fg.id]
  vert.attributes["label"] = lbl #fg.v[fg.id]
  fgl.IDs[lbl] = id # fg.id, to help find it

  # used for cloudgraph solving
  vert.attributes["ready"] = ready
  vert.attributes["backendset"] = 0

  fgl.v[id] = vert # TODO -- this is likely not required, but is used in subgraph methods
  nothing
end

# Add node to graph, given graph struct, labal, init values,
# std dev [TODO -- generalize], particle size and ready flag for concurrency
function addNode!(fg::FactorGraph, lbl::ASCIIString, initval=[0.0]', stdev=[1.0]'; N::Int=100, ready::Int=1)
  fg.id+=1
  vert = ExVertex(fg.id,lbl)
  addNewVarVertInGraph!(fg, vert, fg.id, lbl, ready)
  # dlapi.setupvertgraph!(fg, vert, fg.id, lbl) #fg.v[fg.id]
  dodims = fg.dimID+1
  # TODO -- vert should not loose information here
  setDefaultNodeData!(vert, initval, stdev, dodims, N) #fg.v[fg.id]

  dlapi.addvertex!(fg, vert) #fg.g ##vertr =

  fg.dimID+=size(initval,1) # rows indicate dimensions, move to last dimension
  push!(fg.nodeIDs,fg.id)
  return vert #fg.v[fg.id]
end

# rethink abstraction, maybe closer to CloudGraph use case a better solution
# function addEdge!(g::FGG,n1,n2)
#   edge = dlapi.makeedge(g, n1, n2)
#   dlapi.addedge!(g, edge)
#   # edge = Graphs.make_edge(g, n1, n2)
#   # Graphs.add_edge!(g, edge)
# end

function getVal(vA::Array{Graphs.ExVertex,1})
  len = length(vA)
  vals = Array{Array{Float64,2},1}()
  cols = Array{Int64,1}()
  push!(cols,0)
  rows = Array{Int64,1}()
  for v in vA
      push!(vals, getVal(v))
      c = size(vals[end],2)
      r = size(vals[end],1)
      push!(cols, floor(Int64,c))
      push!(rows, floor(Int64,r))
  end
  cols = cumsum(cols)
  sc = cols[end]
  rw = floor(Int64,rows[1])
  val = Array{Float64,2}(rw, sc)
  for i in 1:(len-1)
      val[:,(cols[i]+1):cols[i+1]] = vals[i]
  end
  val[:,(cols[len]+1):cols[len+1]] = vals[len] # and the last one
  return val
end

# function getDev(v)
#   return v.attributes["stdev"]
# end

function evalFactor(fg::FactorGraph, evalstr::AbstractString)
  return eval(parse(evalstr))
end

function evalFactor(fg::FactorGraph, fct::Graphs.ExVertex)
  return evalFactor(fg, fct.attributes["evalstr"] )
end





function setDefaultFactorNode!(fact::Graphs.ExVertex, f::Union{Pairwise,Singleton})
  # fact.attributes["fnc"] = f
  # fact.attributes["fncargvID"] = Dict{Int,Int}()
  # fact.attributes["eliminated"] = false
  # fact.attributes["potentialused"] = false

  data = FunctionNodeData{typeof(f)}(Int64[], false, false, Int64[], f)
  fact.attributes["data"] = data

  # for graphviz drawing
  fact.attributes["shape"] = "point"  # fg.f[fg.id]
  fact.attributes["width"] = 0.2 #  fg.f[fg.id]

  nothing
end
# function setFncArgIDs!(fact::Graphs.ExVertex, idx::Int64, index::Int64)
#   push!(fact.attributes["data"].fncargvID, index) #[idx]=index
#   #fact.attributes["fncargvID"][idx] = index
#   nothing
# end

function addNewFncVertInGraph!(fgl::FactorGraph, vert::Graphs.ExVertex, id::Int64, lbl::AbstractString, ready::Int)
  vert.attributes = Graphs.AttributeDict() #fg.v[fg.id]
  vert.attributes["label"] = lbl #fg.v[fg.id]
  fgl.f[id] = vert # TODO -- not sure if this is required
  fgl.fIDs[lbl] = id # fg.id

    # used for cloudgraph solving
    vert.attributes["ready"] = ready
    vert.attributes["backendset"] = 0
  nothing
end

function addFactor!(fg::FactorGraph, Xi::Array{Graphs.ExVertex,1},f::Union{Pairwise,Singleton}; ready::Int=1)
  namestring = ""
  for vert in Xi #f.Xi
    namestring = string(namestring,vert.attributes["label"])
  end
  fg.id+=1

  newvert = ExVertex(fg.id,namestring)
  addNewFncVertInGraph!(fg, newvert, fg.id, namestring, ready)
  setDefaultFactorNode!(newvert, f) # fg.f[fg.id]
  push!(fg.factorIDs,fg.id)

  for vert in Xi
    push!(newvert.attributes["data"].fncargvID, vert.index)
  end

  # TODO -- multiple accesses to DB with this method, must refactor!
  newvert = dlapi.addvertex!(fg, newvert)  # used to be two be three lines up ##fg.g
  # idx=0
  for vert in Xi #f.Xi
    dlapi.makeaddedge!(fg, vert, newvert)
    # edge = dlapi.makeedge(fg.g, vert, newvert)
    # dlapi.addedge!(fg.g, edge)
    ## addEdge!(fg.g, vert, newvert) #fg.f[fg.id]
    ## idx+=1
    ## setFncArgIDs!(newvert, idx, vert.index) #fg.f[fg.id]

    #-- push!(newvert.attributes["data"].fncargvID, vert.index)
  end

  return newvert # fg.f[fg.id]
end


function emptyFactorGraph()
    fg = FactorGraph(Graphs.inclist(Graphs.ExVertex,is_directed=false),
                     Graphs.inclist(Graphs.ExVertex,is_directed=true),
                     Dict{Int,Graphs.ExVertex}(),
                     Dict{Int,Graphs.ExVertex}(),
                     Dict{AbstractString,Int}(),
                     Dict{AbstractString,Int}(),
                     0,
                     [],
                     [],
                     Dict{Int,Graphs.ExVertex}(),
                     0,
                     0,
                     nothing,
                     Dict{Int64,Int64}() )
    return fg
end

function prtslperr(s)
  println(s)
  sleep(0.1)
  error(s)
end

# for computing the Bayes Net-----------------------------------------------------
function getEliminationOrder(fg::FactorGraph; ordering::Symbol=:qr)
    s = fg.nodeIDs
    sf = fg.factorIDs
    A=convert(Array{Int},adjacency_matrix(fg.g))[sf,s] # TODO -- order seems brittle
    p = Int[]
    if ordering==:chol
      p = cholfact(A'A,:U,Val{true})[:p] #,pivot=true
    elseif ordering==:qr
      q,r,p = qr(A,Val{true})
    else
      prtslperr("getEliminationOrder -- cannot do the requested ordering $(ordering)")
    end

    # we need the IDs associated with the Graphs.jl and our Type fg
    return fg.nodeIDs[p]
end


# lets create all the vertices first and then deal with the elimination variables thereafter
function addBayesNetVerts!(fg::FactorGraph, elimOrder::Array{Int64,1})
  for p in elimOrder
    vert = localapi.getvertex(fg, p)
    @show vert.label, getData(vert).BayesNetVertID
    if getData(vert).BayesNetVertID == 0   #fg.v[p].
      fg.bnid+=1
      vert.attributes["data"].BayesNetVertID = p # fg.v[p] ##fg.bnverts[p]
      localapi.updatevertex!(fg, vert)
      # fg.bnverts[p] = Graphs.add_vertex!(fg.bn, ExVertex(fg.bnid,string("BayesNet",fg.bnid)))
      # fg.bnverts[p].attributes["label"] = fg.v[p].attributes["label"]
    else
      println("addBayesNetVerts -- something is very wrong, should not have a Bayes net vertex")
    end
  end
end

function addConditional!(fg::FactorGraph, vertID::Int64, lbl, Si)
  bnv = localapi.getvertex(fg, vertID) #fg.v[vertID]
  bnvd = bnv.attributes["data"]
  bnvd.separator = Si
  # bnv.attributes["separator"] = Si # TODO -- remove
  # bnvert = bnv.attributes["BayesNetVert"] # TODO -- remove
  for s in Si
    push!(bnvd.BayesNetOutVertIDs, s)
    # addEdge!(fg.bn, fg.v[s].attributes["BayesNetVert"], bnvert) # TODO -- remove
  end
  localapi.updatevertex!(fg, bnv)
  nothing
end

function addChainRuleMarginal!(fg::FactorGraph, Si)
  lbls = ASCIIString[]
  genmarg = GenericMarginal()
  Xi = Graphs.ExVertex[] #genmarg.Xi = Graphs.ExVertex[]
  for s in Si
    # push!(lbls,fg.v[s].attributes["label"])
    push!(Xi, localapi.getvertex(fg, s))  # push!(Xi,fg.v[s])
  end
  println("adding marginal to")
  for x in Xi @show x.index end
  # addFactor!(fg, marg, lbls)
  addFactor!(fg, Xi, genmarg)
  nothing
end

# TODO -- Cannot have any CloudGraph calls at this level, must refactor
function rmVarFromMarg(fgl::FactorGraph, fromvert::Graphs.ExVertex, gm::Array{Graphs.ExVertex,1})
  for m in gm
    # get all out edges
    # get neighbors
    for n in localapi.outneighbors(fgl, m)
      if n.index == fromvert.index
        alleids = m.attributes["data"].edgeIDs
        # println("consider marginal, edge of interest between $(m.index) and $(n.index)")
        i = 0
        for id in alleids
          # println("rmVarFromMarg -- at id=$(id)")
          i+=1
          edge = localapi.getedge(fgl, id)
          if edge != nothing # hack to avoid dictionary use case
            if edge.SourceVertex.exVertexId == m.index || edge.DestVertex.exVertexId == m.index
              warn("removing edge $(edge.neo4jEdgeId), between $(m.index) and $(n.index)")
              localapi.deleteedge!(fgl, edge)
              m.attributes["data"].edgeIDs = alleids[[collect(1:(i-1));collect((i+1):length(alleids))]]
              localapi.updatevertex!(fgl, m)
            end
          end
        end
      end
    end
    # if 0 edges, delete the marginal
    if length(localapi.outneighbors(fgl, m)) <= 1
      warn("removing vertex id=$(m.index)")
      localapi.deletevertex!(fgl,m)
    end
  end
  nothing
end

function buildBayesNet!(fg::FactorGraph, p::Array{Int,1})
    addBayesNetVerts!(fg, p)
    for v in p
      println()
      println("Eliminating $(v)")
      println("===============")
      println()
      # which variable are we eliminating
      #fg.v[v].attributes["label"]

      # all factors adjacent to this variable
      fi = Int64[]
      Si = Int64[]
      gm = ExVertex[]
      # TODO -- optimize outneighbor calls like this
      vert = localapi.getvertex(fg,v)
      for fct in localapi.outneighbors(fg, vert) #out_neighbors(dlapi.getvertex(fg,v),fg.g) # (fg.v[v]
        if (getData(fct).eliminated != true) # && fct.attributes["ready"]==1 && fct.attributes["backendset"]==1 #if (fct.attributes["eliminated"] != true)
          push!(fi, fct.index)
          for sepNode in localapi.outneighbors(fg, fct) #out_neighbors(fct,fg.g)
            if sepNode.index != v && length(findin(sepNode.index,Si)) == 0
              push!(Si,sepNode.index)
            end
          end
          fct.attributes["data"].eliminated = true #fct.attributes["eliminated"] = true
          localapi.updatevertex!(fg, fct) # TODO -- this might be a premature statement
        end

        if typeof(getData(fct).fnc) == GenericMarginal
          push!(gm, fct)
        end
      end

      if v != p[end]
        addConditional!(fg, v, "", Si)
        # not yet inserting the new prior p(Si) back into the factor graph
      end

      tuv = localapi.getvertex(fg, v) # TODO -- This may well through away valuable data
      tuv.attributes["data"].eliminated = true # fg.v[v].
      localapi.updatevertex!(fg, tuv)
      # dlapi.getvertex(fg,v).attributes["data"].eliminated = true # fg.v[v].
      ## fg.v[v].attributes["eliminated"] = true


      # TODO -- remove links from current vertex to any marginals
      rmVarFromMarg(fg, vert, gm)

      #add marginal on remaining variables... ? f(xyz) = f(x | yz) f(yz)
      # new function between all Si
      addChainRuleMarginal!(fg, Si)

      # readline(STDIN)
    end
    nothing
end

# some plotting functions on the factor graph
function stackVertXY(fg::FactorGraph, lbl::ASCIIString)
    v = dlapi.getvertex(fg,lbl) # v = fg.v[fg.IDs[lbl]]
    vals = getVal(v)
    X=vec(vals[1,:])
    Y=vec(vals[2,:])
    # X=vec(v.attributes["val"][1,:])
    # Y=vec(v.attributes["val"][2,:])
    return X,Y
end

function getKDE(v::Graphs.ExVertex)
  #return kde!(getVal(v), "lcv") # TODO - getBW(v)
  return kde!(getVal(v), getBWVal(v)[:,1])
end

function getVertKDE(fgl::FactorGraph, id::Int64)
  v = dlapi.getvertex(fgl,id)  # v = fgl.v[fgl.IDs[lbl]]
  return getKDE(v)
end
function getVertKDE(fgl::FactorGraph, lbl::ASCIIString)
  v = dlapi.getvertex(fgl,lbl)  # v = fgl.v[fgl.IDs[lbl]]
  return getKDE(v)
end


function drawCopyFG(fgl::FactorGraph)
  fgd = deepcopy(fgl)
  for v in fgd.v
    delete!(v[2].attributes,"data")
    # delete!(v[2].attributes,"initstdev")
    # delete!(v[2].attributes,"dimIDs")
    # delete!(v[2].attributes,"eliminated")
    # delete!(v[2].attributes,"val")
    # delete!(v[2].attributes,"bw")
    # delete!(v[2].attributes,"BayesNetVert")
    # delete!(v[2].attributes,"initval")
    # delete!(v[2].attributes,"dims")
  end
  for v in fgd.f
    delete!(v[2].attributes,"data")
    # delete!(v[2].attributes,"fncargvID")
    # delete!(v[2].attributes,"eliminated")
    # delete!(v[2].attributes,"fnc")
    # delete!(v[2].attributes,"potentialused")
  end
  return fgd
end

function writeGraphPdf(fgl::FactorGraph)
  fgd = drawCopyFG(fgl)
  println("Writing factor graph file")
  fid = open("fg.dot","w+")
  write(fid,to_dot(fgd.g))
  close(fid)
  run(`dot fg.dot -Tpdf -o fg.pdf`)
  nothing
end




function appendUseFcts!(usefcts, lblid::Int, fct::Graphs.ExVertex, fid::Int)
  for tp in usefcts
    if tp[2].label == fct.label
      # println("Skipping repeat of $(fct.label)")
      return nothing
    end
  end
  tpl = (fct.index, fct) #(lblid, fct, fid)
  push!(usefcts, tpl )
  nothing
end


function expandEdgeListNeigh!(fgl::FactorGraph,
                              vertdict::Dict{Int64,Graphs.ExVertex},
                              edgedict::Dict{Int64,Graphs.Edge{Graphs.ExVertex}})
  #asfd
  for vert in vertdict
    for newedge in out_edges(vert[2],fgl.g)
      if !haskey(edgedict, newedge.index)
        edgedict[newedge.index] = newedge
      end
    end
  end

  nothing
end

# dictionary of unique vertices from edgelist
function expandVertexList!(fgl::FactorGraph,
  edgedict::Dict{Int64,Graphs.Edge{Graphs.ExVertex}},
  vertdict::Dict{Int64,Graphs.ExVertex})

  # go through all source and target nodes
  for edge in edgedict
    if !haskey(vertdict, edge[2].source.index)
      vertdict[edge[2].source.index] = edge[2].source
    end
    if !haskey(vertdict, edge[2].target.index)
      vertdict[edge[2].target.index] = edge[2].target
    end
  end
  nothing
end

function edgelist2edgedict(edgelist::Array{Graphs.Edge{Graphs.ExVertex},1})
  edgedict = Dict{Int64,Graphs.Edge{Graphs.ExVertex}}()
  for edge in edgelist
    edgedict[edge.index] = edge
  end
  return edgedict
end

function addVerticesSubgraph(fgl::FactorGraph,
    fgseg::FactorGraph,
    vertdict::Dict{Int64,Graphs.ExVertex})

    for vert in vertdict
      fgseg.g.vertices[vert[1]] = vert[2]
      if haskey(fgl.v,vert[1])
        fgseg.v[vert[1]] = vert[2]
        fgseg.IDs[vert[2].label] = vert[1]

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
        fgseg.fIDs[vert[2].label] = vert[1]
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

# NOTICE, nodeIDs and factorIDs are not brough over by this method yet
# must sort out for incremental updates
function genSubgraph(fgl::FactorGraph,
    vertdict::Dict{Int64,Graphs.ExVertex})
    # edgedict::Dict{Int64,Graphs.Edge{Graphs.ExVertex}},

  fgseg = FactorGraph() # new handle for just a segment of the graph
  fgseg.g = Graphs.inclist(Graphs.ExVertex,is_directed=false)

  fgseg.v = Dict{Int,Graphs.ExVertex}()
  fgseg.f = Dict{Int,Graphs.ExVertex}()
  fgseg.IDs = Dict{AbstractString,Int}()
  fgseg.fIDs = Dict{AbstractString,Int}()

  fgseg.g.vertices = Array{Graphs.ExVertex,1}(length(fgl.g.vertices))
  fgseg.g.inclist = Array{Array{Graphs.Edge{Graphs.ExVertex},1},1}(length(fgl.g.inclist))

  addVerticesSubgraph(fgl, fgseg, vertdict)

  fgseg.id = fgl.id
  fgseg.bnid = fgl.bnid
  fgseg.dimID = fgl.dimID

  return fgseg
end

function getShortestPathNeighbors(fgl::FactorGraph;
    from::Graphs.ExVertex=nothing,
    to::Graphs.ExVertex=nothing,
    neighbors::Int=0 )

  edgelist = shortest_path(fgl.g, ones(num_edges(fgl.g)), from, to)
  vertdict = Dict{Int64,Graphs.ExVertex}()
  edgedict = edgelist2edgedict(edgelist)
  expandVertexList!(fgl, edgedict, vertdict) # grow verts
  for i in 1:neighbors
    # println("Expanding subgraph edge dict")
    expandEdgeListNeigh!(fgl, vertdict, edgedict) # grow edges
    expandVertexList!(fgl, edgedict, vertdict) # grow verts
  end
  return vertdict
end

function subgraphShortestPath(fgl::FactorGraph;
    from::Graphs.ExVertex=nothing,
    to::Graphs.ExVertex=nothing,
    neighbors::Int=0  )

  vertdict = getShortestPathNeighbors(fgl, from=from, to=to, neighbors=neighbors)
  return genSubgraph(fgl, vertdict)
end

# explore all shortest paths combinations in verts, add neighbors and reference subgraph
function subgraphFromVerts(fgl::FactorGraph,
    verts::Dict{Int64,Graphs.ExVertex};
    neighbors::Int=0  )

    allverts = Dict{Int64,Graphs.ExVertex}()
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
    verts::Array{Int64,1};
    neighbors::Int=0  )

  vertdict = Dict{Int,Graphs.ExVertex}()
  for vert in verts
    vertdict[vert] = fgl.v[vert]
  end

  return subgraphFromVerts(fgl,vertdict,neighbors=neighbors)
end

# explore all shortest paths combinations in verts, add neighbors and reference subgraph
# Using unique index into graph data structure
function subgraphFromVerts(fgl::FactorGraph,
    verts::Array{ASCIIString,1};
    neighbors::Int=0  )

  vertdict = Dict{Int,Graphs.ExVertex}()
  for vert in verts
    id = -1
    if haskey(fgl.IDs, vert)
      id = fgl.IDs[vert]
    else
      error("FactorGraph01 only supports variable node subgraph search at this point")
    end
    vertdict[id] = fgl.v[id]
  end

  return subgraphFromVerts(fgl,vertdict,neighbors=neighbors)
end
