import Base.convert
import Base.==

abstract Pairwise
abstract Singleton

reshapeVec2Mat(vec::Vector, rows::Int) = reshape(vec, rows, round(Int,length(vec)/rows))
# function reshapeVec2Mat(vec::Vector, rows::Int)
#   M = reshape(vec, rows, round(Int,length(vec)/rows))
#   return ndims(M) < 2 ? (M')' : M
# end

type FactorGraph
  g
  bn
  v::Dict{Int,Graphs.ExVertex}
  f::Dict{Int,Graphs.ExVertex}
  IDs::Dict{AbstractString,Int}
  fIDs::Dict{AbstractString,Int}
  id::Int
  nodeIDs::Array{Int,1}
  factorIDs::Array{Int,1}
  bnverts::Dict{Int,Graphs.ExVertex}
  bnid::Int
  dimID::Int64
end

type VariableNodeData
  initval::Array{Float64,2}
  initstdev::Array{Float64,2}
  val::Array{Float64,2}
  bw::Array{Float64,2}
  BayesNetOutVertIDs::Array{Int64,1}
  dimIDs::Array{Int64,1}
  dims::Int64
  eliminated::Bool
  BayesNetVertID::Int64
  separator::Array{Int64,1} # Will bring to hard type soon, don't worry about this one just yet
  VariableNodeData() = new()
  VariableNodeData(x...) = new(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10])
end

type PackedVariableNodeData
  vecinitval::Array{Float64,1}
  diminitval::Int64
  vecinitstdev::Array{Float64,1}
  diminitdev::Int64
  vecval::Array{Float64,1}
  dimval::Int64
  vecbw::Array{Float64,1}
  dimbw::Int64
  BayesNetOutVertIDs::Array{Int64,1}
  dimIDs::Array{Int64,1}
  dims::Int64
  eliminated::Bool
  BayesNetVertID::Int64
  separator::Array{Int64,1}
  PackedVariableNodeData() = new()
  PackedVariableNodeData(x...) = new(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14])
end


function convert(::Type{PackedVariableNodeData}, d::VariableNodeData)
  return PackedVariableNodeData(d.initval[:],size(d.initval,1),
                              d.initstdev[:],size(d.initstdev,1),
                              d.val[:],size(d.val,1),
                              d.bw[:],size(d.bw,1),
                              d.BayesNetOutVertIDs,
                              d.dimIDs, d.dims, d.eliminated,
                              d.BayesNetVertID, d.separator)
end
function convert(::Type{VariableNodeData}, d::PackedVariableNodeData)
  r1 = d.diminitval
  c1 = floor(Int,length(d.vecinitval)/r1)
  M1 = reshape(d.vecinitval,r1,c1)

  r2 = d.diminitdev
  c2 = floor(Int,length(d.vecinitstdev)/r2)
  M2 = reshape(d.vecinitstdev,r2,c2)

  r3 = d.dimval
  c3 = floor(Int,length(d.vecval)/r3)
  M3 = reshape(d.vecval,r3,c3)

  r4 = d.dimbw
  c4 = floor(Int,length(d.vecbw)/r4)
  M4 = reshape(d.vecbw,r4,c4)

  return VariableNodeData(M1,M2,M3,M4, d.BayesNetOutVertIDs,
    d.dimIDs, d.dims, d.eliminated, d.BayesNetVertID, d.separator)
end

function compare(a::VariableNodeData,b::VariableNodeData)
    TP = true
    TP = TP && a.initval == b.initval
    TP = TP && a.initstdev == b.initstdev
    TP = TP && a.val == b.val
    TP = TP && a.bw == b.bw
    TP = TP && a.BayesNetOutVertIDs == b.BayesNetOutVertIDs
    TP = TP && a.dimIDs == b.dimIDs
    TP = TP && a.dims == b.dims
    TP = TP && a.eliminated == b.eliminated
    TP = TP && a.BayesNetVertID == b.BayesNetVertID
    TP = TP && a.separator == b.separator
    return TP
end

function ==(a::VariableNodeData,b::VariableNodeData)
  return compare(a,b)
end



function getVal(v::Graphs.ExVertex)
  return v.attributes["data"].val
  # return v.attributes["val"]
end
function setVal!(v::Graphs.ExVertex, val::Array{Float64,2})
  v.attributes["data"].val = val
  # v.attributes["val"] = val
  nothing
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

# function setDefaultNodeDataOld!(v::Graphs.ExVertex, initval::Array{Float64,2},
#                               stdev::Array{Float64,2}, dodims::Int64, N::Int64)
#   pN = Union{}
#   if size(initval,2) < N
#     p = kde!(initval,diag(stdev));
#     pN = resample(p,N)
#   else
#     pN = kde!(initval, "lcv")
#   end
#   # v.attributes["initval"] = initval
#   # v.attributes["initstdev"] = stdev
#   # v.attributes["eliminated"] = false
#   v.attributes["BayesNetVert"] = Union{}
#   dims = size(initval,1) # rows indicate dimensions
#   # v.attributes["dims"] = dims
#   # v.attributes["dimIDs"] = round(Int64,linspace(dodims,dodims+dims-1,dims))
#
#   setDefaultNodeDataReplace!(v,initval,stdev,dodims,N)
#   # setVal!(v, getPoints(pN), getBW(pN)[:,1])
#   nothing
# end

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
  # v.attributes["BayesNetVert"] = Union{} # TODO -- remove
  nothing
end

function addNode!(fg::FactorGraph, lbl, initval=[0.0]', stdev=[1.0]'; N::Int=100)
  fg.id+=1
  fg.v[fg.id] = dlapi.addvertex!(fg.g, ExVertex(fg.id,lbl))
  # fg.v[fg.id] = Graphs.add_vertex!(fg.g, ExVertex(fg.id,lbl))
  fg.IDs[lbl] = fg.id
  fg.v[fg.id].attributes = Graphs.AttributeDict()
  fg.v[fg.id].attributes["label"] = lbl
  dodims = fg.dimID+1

  setDefaultNodeData!(fg.v[fg.id], initval, stdev, dodims, N)

  fg.dimID+=size(initval,1) # rows indicate dimensions, move to last dimension
  push!(fg.nodeIDs,fg.id)
  return fg.v[fg.id]
end

function addEdge!(g,n1,n2)
  edge = dlapi.makeedge(g, n1, n2)
  dlapi.addedge!(g, edge)
  # edge = Graphs.make_edge(g, n1, n2)
  # Graphs.add_edge!(g, edge)
end

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

# function FactorEvalStr(fg,factor)
#   evalstr = string(factor.attributes["data"].fnc) #string(factor.attributes["fnc"])
#   evalstr = string(evalstr,"(")
#   dd = factor.attributes["data"].fncargvID #attributes["fncargvID"]
#   idx = 0
#   for nodeid in dd
#     idx += 1
#     evalstr=string(evalstr,"getVal(fg.v[",nodeid,"]),")
#     # evalstr=string(evalstr,"getVal(fg.v[",dd[idx],"]),")
#   end
#   evalstr = string(chop(evalstr),")")
#   return evalstr
# end

function evalFactor(fg::FactorGraph, evalstr::AbstractString)
  return eval(parse(evalstr))
end

function evalFactor(fg::FactorGraph, fct::Graphs.ExVertex)
  return evalFactor(fg, fct.attributes["evalstr"] )
end


type FunctionNodeData{T}
  fncargvID::Array{Int64,1}
  eliminated::Bool
  potentialused::Bool
  fnc::T
end


# function ohdear(customtype)
#   vertex.packed
#   vertex.templatetype = "MyType"
#   SPT = parse..(vertextemplatetype)
#   readproto(convert(Vector{UInt8},props["packed"]), customtype{SPT}() ) # something to this effect
# end

function setDefaultFactorNode!(fact::Graphs.ExVertex, f::Union{Pairwise,Singleton})
  # fact.attributes["fnc"] = f
  # fact.attributes["fncargvID"] = Dict{Int,Int}()
  # fact.attributes["eliminated"] = false
  # fact.attributes["potentialused"] = false

  data = FunctionNodeData{typeof(f)}(Int64[], false, false, f)
  fact.attributes["data"] = data
  nothing
end
function setFncArgIDs!(fact::Graphs.ExVertex, idx::Int64, index::Int64)
  push!(fact.attributes["data"].fncargvID, index) #[idx]=index
  #fact.attributes["fncargvID"][idx] = index
  nothing
end

function addFactor!(fg::FactorGraph, Xi::Array{Graphs.ExVertex,1},f::Union{Pairwise,Singleton})
  namestring = ""
  for vert in Xi #f.Xi
    namestring = string(namestring,vert.attributes["label"])
  end
  fg.id+=1
  fg.f[fg.id] = dlapi.addvertex!(fg.g, ExVertex(fg.id,namestring))
  # fg.f[fg.id] = Graphs.add_vertex!(fg.g, ExVertex(fg.id,namestring))
  fg.fIDs[namestring] = fg.id
  # for graphviz drawing
  fg.f[fg.id].attributes["shape"] = "point"
  fg.f[fg.id].attributes["width"] = 0.2

  setDefaultFactorNode!(fg.f[fg.id], f)
  push!(fg.factorIDs,fg.id)

  idx=0
  for vert in Xi #f.Xi
    addEdge!(fg.g,vert,fg.f[fg.id])
    idx+=1
    setFncArgIDs!(fg.f[fg.id],idx,vert.index)
  end

  #fg.f[fg.id].attributes["evalstr"] = FactorEvalStr(fg,fg.f[fg.id])
  return fg.f[fg.id]
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
                     0)
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
    A=convert(Array{Int},adjacency_matrix(fg.g))[sf,s]
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
# function addBayesNetVertsOld!(fg::FactorGraph, elimOrder::Array{Int64,1})
#   for p in elimOrder
#     if (fg.v[p].attributes["BayesNetVert"] == Union{})
#       fg.bnid+=1
#       fg.bnverts[p] = dlapi.addvertex!(fg.bn, ExVertex(fg.bnid,string("BayesNet",fg.bnid)))
#       fg.v[p].attributes["BayesNetVert"] = fg.bnverts[p]
#       fg.bnverts[p].attributes["label"] = fg.v[p].attributes["label"]
#     else
#       println("addBayesNetVerts -- something is very wrong, should not have a Bayes net vertex")
#     end
#   end
# end

# lets create all the vertices first and then deal with the elimination variables thereafter
function addBayesNetVerts!(fg::FactorGraph, elimOrder::Array{Int64,1})
  for p in elimOrder
    if fg.v[p].attributes["data"].BayesNetVertID == 0
      fg.bnid+=1
      # fg.bnverts[p] = Graphs.add_vertex!(fg.bn, ExVertex(fg.bnid,string("BayesNet",fg.bnid)))
      fg.v[p].attributes["data"].BayesNetVertID = p#fg.bnverts[p]
      # fg.bnverts[p].attributes["label"] = fg.v[p].attributes["label"]
    else
      println("addBayesNetVerts -- something is very wrong, should not have a Bayes net vertex")
    end
  end
end

function addConditional!(fg::FactorGraph, vertID, lbl, Si)
  bnv = fg.v[vertID]
  bnvd = bnv.attributes["data"]
  bnvd.separator = Si
  # bnv.attributes["separator"] = Si # TODO -- remove
  # bnvert = bnv.attributes["BayesNetVert"] # TODO -- remove
  for s in Si
    push!(bnvd.BayesNetOutVertIDs, s)
    # addEdge!(fg.bn, fg.v[s].attributes["BayesNetVert"], bnvert) # TODO -- remove
  end
end

function addChainRuleMarginal!(fg::FactorGraph, Si)
  lbls = ASCIIString[]
  genmarg = GenericMarginal()
  Xi = Graphs.ExVertex[] #genmarg.Xi = Graphs.ExVertex[]
  for s in Si
    # push!(lbls,fg.v[s].attributes["label"])
    push!(Xi,fg.v[s])
  end
  # addFactor!(fg, marg, lbls)
  addFactor!(fg,Xi, genmarg)
  Union{}
end

function buildBayesNet!(fg::FactorGraph, p::Array{Int,1})
    addBayesNetVerts!(fg, p)
    for v in p
      # which variable are we eliminating
      #fg.v[v].attributes["label"]

      # all factors adjacent to this variable
      fi = Int64[]
      Si = Int64[]
      for fct in out_neighbors(fg.v[v],fg.g)
        if (fct.attributes["data"].eliminated != true) #if (fct.attributes["eliminated"] != true)
          push!(fi, fct.index)
          for sepNode in out_neighbors(fct,fg.g)
            if (sepNode.index != v && length(findin(sepNode.index,Si)) == 0)
              push!(Si,sepNode.index)
            end
          end
          fct.attributes["data"].eliminated = true #fct.attributes["eliminated"] = true
        end
      end
      #@show fg.v[v].attributes["sepfactors"] = fi
      addConditional!(fg, v, "", Si)
      # not yet inserting the new prior p(Si) back into the factor graph
      fg.v[v].attributes["data"].eliminated = true
      # fg.v[v].attributes["eliminated"] = true

      #add marginal on remaining variables... ? f(xyz) = f(x | yz) f(yz)
      # new function between all Si
      addChainRuleMarginal!(fg, Si)
    end
    nothing
end

# BayesTree declarations
type BayesTree
  bt
  btid::Int
  cliques::Dict{Int,Graphs.ExVertex}
  #edges
  frontals::Dict{ASCIIString,Int64}
end

function emptyBayesTree()
    bt =   BayesTree(Graphs.inclist(Graphs.ExVertex,is_directed=true),
                     0,
                     Dict{Int,Graphs.ExVertex}(),
                     #[],
                     Dict{AbstractString, Int}())
    return bt
end

# create a new clique
function addClique!(bt::BayesTree, fg::FactorGraph, varID::Int, condIDs::Array{Int}=Int[])
  bt.btid += 1
  clq = Graphs.add_vertex!(bt.bt, ExVertex(bt.btid,string("Clique",bt.btid)))
  bt.cliques[bt.btid] = clq
  clq.attributes["frontalIDs"] = Int[]
  clq.attributes["conditIDs"] = Int[]

  clq.attributes["label"] = ""
  clq.attributes["potentials"] = []

  appendClique!(bt, bt.btid, fg, varID, condIDs)
  return clq
end

# generate the label for particular clique -- graphviz drawing
function makeCliqueLabel(fgl::FactorGraph, bt::BayesTree, clqID::Int)
  clq = bt.cliques[clqID]
  flbl = ""
  clbl = ""
  for fr in clq.attributes["frontalIDs"]
    flbl = string(flbl,fgl.v[fr].attributes["label"], ",")
  end
  for cond in clq.attributes["conditIDs"]
    clbl = string(clbl,fgl.v[cond].attributes["label"], ",")
  end
  clq.attributes["label"] = string(flbl, ": ", clbl)
end

# add a conditional ID to clique
function appendConditional(bt::BayesTree, clqID::Int, condIDs::Array{Int,1})
  #@show "adding conditionals", condIDs
  clq = bt.cliques[clqID]
  clq.attributes["conditIDs"] = union(clq.attributes["conditIDs"], condIDs)
end

# Add a new frontal variable to clique
function appendClique!(bt::BayesTree, clqID::Int, fg::FactorGraph, varID::Int, condIDs::Array{Int,1}=Int[])
  clq = bt.cliques[clqID]
  var = fg.v[varID]
  # add frontal variable
  push!(clq.attributes["frontalIDs"],varID)
  # total dictionary of frontals for easy access
  bt.frontals[var.attributes["label"]] = clqID#bt.btid

  appendConditional(bt, clqID, condIDs)
  makeCliqueLabel(fg, bt, clqID)
  Union{}
end


# instantiate a new child clique in the tree
function newChildClique!(bt::BayesTree, fg::FactorGraph, CpID::Int, varID::Int, Sepj::Array{Int,1})
  chclq = addClique!(bt, fg, varID, Sepj)
  parent = bt.cliques[CpID]
  addEdge!(bt.bt, parent, chclq)

  return chclq
end

# post order tree traversal and build potential functions
function findCliqueFromFrontal(bt::BayesTree, frtlID::Int)
  for cliqPair in bt.cliques
    id = cliqPair[1]
    cliq = cliqPair[2]
    for frtl in cliq.attributes["frontalIDs"]
      if frtl == frtlID
        return cliq
      end
    end
  end
  error("Clique with desired frontal ID not found")
end


# eliminate a variable for new
function newPotential(tree::BayesTree, fg::FactorGraph, var::Int, prevVar::Int, p::Array{Int,1})
    #@show fg.v[var].attributes["label"]
    if (length(fg.v[var].attributes["data"].separator) == 0) #(length(fg.v[var].attributes["separator"]) == 0)
      if (length(tree.cliques) == 0)
        addClique!(tree, fg, var)
      else
        appendClique!(tree, 1, fg, var) # add to root
      end
    else
      Sj = fg.v[var].attributes["data"].separator
      # Sj = fg.v[var].attributes["separator"]
      # find parent clique Cp that containts the first eliminated variable of Sj as frontal
      firstelim = 99999999999
      for s in Sj
        temp = findfirst(p, s)
        if (temp < firstelim)
          firstelim = temp
        end
      end
      felbl = fg.v[p[firstelim]].attributes["label"]
      CpID = tree.frontals[felbl]
      # look to add this conditional to the tree
      unFC = union(tree.cliques[CpID].attributes["frontalIDs"], tree.cliques[CpID].attributes["conditIDs"])
      if (sort(unFC) == sort(Sj))
        appendClique!(tree, CpID, fg, var)
      else
        newChildClique!(tree, fg, CpID, var, Sj)
      end
    end
end

# build the whole tree in batch format
function buildTree!(tree::BayesTree, fg::FactorGraph, p::Array{Int,1})
  rp = flipdim(p,1)
  prevVar = 0
  for var in rp
    newPotential(tree, fg, var, prevVar, p)
    prevVar = var
  end
end


marg(x...) = -1, [0.0] # TODO -- this should be removed


## Find batch belief propagation solution
function prepBatchTree!(fg::FactorGraph; ordering::Symbol=:qr,drawpdf::Bool=false)
  p = getEliminationOrder(fg, ordering=ordering)
  #@show p=[1, 3, 8, 10, 5]
  # for v in p
  #     print("$(fg.v[v].label),")
  # end
  println()

  fge = deepcopy(fg)
  buildBayesNet!(fge, p)

  tree = emptyBayesTree()
  buildTree!(tree, fge, p)

  println("Bayes Net")
  sleep(0.1)
  #fid = open("bn.dot","w+")
  #write(fid,to_dot(fge.bn))
  #close(fid)

  #GraphViz.Graph(to_dot(fge.bn))
  # Michael reference -- x2->x1, x2->x3, x2->x4, x2->l1, x4->x3, l1->x3, l1->x4
  println("Bayes Tree")
  if drawpdf
    fid = open("bt.dot","w+")
    write(fid,to_dot(tree.bt))
    close(fid)
    run(`dot bt.dot -Tpdf -o bt.pdf`)
  end

  # GraphViz.Graph(to_dot(tree.bt))
  #Michael reference 3sig -- x2l1x4x3    x1|x2

  println("Find potential functions for each clique")
  cliq = tree.cliques[1] # start at the root
  buildCliquePotentials(fg, tree, cliq);

  return tree
end

function resetFactorGraphNewTree!(fg::FactorGraph)
  for v in fg.v
    v[2].attributes["data"].eliminated = false
    # v[2].attributes["eliminated"] = false
    v[2].attributes["data"].BayesNetOutVertIDs = Int64[]
    v[2].attributes["data"].BayesNetVertID = 0
    v[2].attributes["data"].separator = Int64[]
    # v[2].attributes["BayesNetVert"] = Union{}
    # v[2].attributes["separator"] = Int64[]
  end
  for f in fg.f
    f[2].attributes["data"].eliminated = false #f[2].attributes["eliminated"] = false
    f[2].attributes["data"].potentialused = false #f[2].attributes["potentialused"] = false
  end

  nothing
end

function wipeBuildNewTree!(fg::FactorGraph; ordering=:qr)
  resetFactorGraphNewTree!(fg);
  return prepBatchTree!(fg, ordering=ordering);
end

function whichCliq(bt::BayesTree, frt::ASCIIString)
    bt.cliques[bt.frontals[frt]]
end

# some plotting functions on the factor graph
function stackVertXY(fg::FactorGraph, lbl::ASCIIString)
    v = fg.v[fg.IDs[lbl]]
    vals = getVal(v)
    X=vec(vals[1,:])
    Y=vec(vals[2,:])
    # X=vec(v.attributes["val"][1,:])
    # Y=vec(v.attributes["val"][2,:])
    return X,Y
end

function getKDE(v::Graphs.ExVertex)
  return kde!(getVal(v), "lcv") # TODO - getBW(v)
end

function getVertKDE(fgl::FactorGraph, lbl::ASCIIString)
  v = fgl.v[fgl.IDs[lbl]]
  return getKDE(v)
end

#
# function get2DSamples(fg::FactorGraph, sym; from::Int64=0, to::Int64=999999999, minnei::Int64=0)
#   X = Array{Float64,1}()
#   Y = Array{Float64,1}()
#
#   # if sym = 'l', ignore single measurement landmarks
#
#   for id in fg.IDs
#       if id[1][1] == sym
#         val = parse(Int,id[1][2:end])
#         if from <= val && val <= to
#           if length(out_neighbors(fg.v[id[2]],fg.g)) >= minnei
#               X=[X;vec(fg.v[id[2]].attributes["val"][1,:])]
#               Y=[Y;vec(fg.v[id[2]].attributes["val"][2,:])]
#           end
#         end
#       end
#   end
#   return X,Y
# end
#
# function getAll2D(fg, sym; minnei::Int64=0)
#   return get2DSamples(fg, sym, minnei=minnei)
# end

#
# function get2DSampleMeans(fg::FactorGraph, sym; from::Int64=0, to::Int64=9999999999, minnei::Int64=0)
#   X = Array{Float64,1}()
#   Y = Array{Float64,1}()
#   Th = Array{Float64,1}()
#   LB = ASCIIString[]
#
#   # if sym = 'l', ignore single measurement landmarks
#   allIDs = Array{Int,1}()
#   for id in fg.IDs
#       if id[1][1] == sym
#           val = parse(Int,id[1][2:end])
#           if from <= val && val <= to
#             if length(out_neighbors(fg.v[id[2]],fg.g)) >= minnei
#               push!(allIDs, val)
#             end
#           end
#       end
#   end
#   allIDs = sort(allIDs)
#
#   for id in allIDs
#       X=[X;Base.mean(vec(fg.v[fg.IDs[string(sym,id)]].attributes["val"][1,:]))]
#       Y=[Y;Base.mean(vec(fg.v[fg.IDs[string(sym,id)]].attributes["val"][2,:]))]
#       if sym == 'x'
#         Th=[Th;Base.mean(vec(fg.v[fg.IDs[string(sym,id)]].attributes["val"][3,:]))]
#       end
#       push!(LB, string(sym,id))
#   end
#   return X,Y,Th,LB
# end
#
# #draw landmark positions
# function getAll2DMeans(fg, sym)
#   return get2DSampleMeans(fg, sym)
# end
#
# function getAll2DPoses(fg::FactorGraph)
#     return getAll2D(fg, 'x')
# end
#
# function get2DPoseSamples(fg::FactorGraph; from::Int64=0, to::Int64=999999999)
#   return get2DSamples(fg, 'x'; from=from, to=to)
# end
#
# function get2DPoseMeans(fg::FactorGraph; from::Int64=0, to::Int64=999999999)
#   return get2DSampleMeans(fg, 'x', from=from, to=to)
# end
#
#
# function get2DPoseMax(fgl::FactorGraph;
#             from::Int=-99999999999, to::Int=9999999999 )
#   xLB,ll = ls(fgl) # TODO add: from, to, special option 'x'
#   X = Array{Float64,1}()
#   Y = Array{Float64,1}()
#   Th = Array{Float64,1}()
#   LB = ASCIIString[]
#   for lbl in xLB
#     if from <= parse(Int,lbl[2:end]) <=to
#       mv = getKDEMax(getVertKDE(fgl,lbl))
#       push!(X,mv[1])
#       push!(Y,mv[2])
#       push!(Th,mv[3])
#       push!(LB, lbl)
#     end
#   end
#   return X, Y, Th, LB
# end
#
# function getAll2DLandmarks(fg::FactorGraph, minnei::Int64=0)
#     return getAll2D(fg, 'l', minnei=minnei)
# end
#
# function get2DLandmSamples(fg::FactorGraph; from::Int64=0, to::Int64=999999999, minnei::Int64=0)
#   return get2DSamples(fg, 'l', from=from, to=to, minnei=minnei)
# end
#
# function get2DLandmMeans(fg::FactorGraph; from::Int64=0, to::Int64=999999999, minnei::Int64=0)
#   return get2DSampleMeans(fg, 'l', from=from, to=to, minnei=minnei)
# end
#
# function removeKeysFromArr(fgl::FactorGraph, torm::Array{Int,1}, lbl::Array{ASCIIString,1})
#   retlbs = ASCIIString[]
#   for i in 1:length(lbl)
#     id = parse(Int,lbl[i][2:end])
#     if findfirst(torm,id) == 0
#       push!(retlbs, lbl[i])
#     else
#       println("removeKeysFromArr -- skipping $(lbl[i]), id=$(id)")
#     end
#   end
#   return retlbs
# end
#
# function get2DLandmMax(fgl::FactorGraph;
#                 from::Int=-99999999999, to::Int=9999999999, showmm=false,MM=Union{} )
#   xLB,lLB = ls(fgl) # TODO add: from, to, special option 'x'
#   if !showmm lLB = removeKeysFromArr(fgl, collect(keys(MM)), lLB); end
#   X = Array{Float64,1}()
#   Y = Array{Float64,1}()
#   Th = Array{Float64,1}()
#   LB = ASCIIString[]
#   for lbl in lLB
#     if from <= parse(Int,lbl[2:end]) <=to
#       mv = getKDEMax(getVertKDE(fgl,lbl))
#       push!(X,mv[1])
#       push!(Y,mv[2])
#       push!(LB, lbl)
#     end
#   end
#   return X, Y, Th, LB
# end

function drawCopyFG(fgl::FactorGraph)
  fgd = deepcopy(fgl)
  for v in fgd.v
    delete!(v[2].attributes,"initstdev")
    delete!(v[2].attributes,"dimIDs")
    delete!(v[2].attributes,"eliminated")
    delete!(v[2].attributes,"val")
    delete!(v[2].attributes,"bw")
    delete!(v[2].attributes,"BayesNetVert")
    delete!(v[2].attributes,"initval")
    delete!(v[2].attributes,"dims")
  end
  for v in fgd.f
    delete!(v[2].attributes,"fncargvID")
    delete!(v[2].attributes,"eliminated")
    delete!(v[2].attributes,"fnc")
    delete!(v[2].attributes,"potentialused")
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
      println("Skipping repeat of $(fct.label)")
      return nothing
    end
  end
  tpl = (fct.index, fct) #(lblid, fct, fid)
  push!(usefcts, tpl )
  nothing
end

function getCliquePotentials!(fg::FactorGraph, bt::BayesTree, cliq::Graphs.ExVertex)
    frtl = cliq.attributes["frontalIDs"]
    cond = cliq.attributes["conditIDs"]
    allids = [frtl;cond]
    alldimIDs = Int[]
    for fid in frtl
      alldimIDs = [alldimIDs;fg.v[fid].attributes["data"].dimIDs]
      # alldimIDs = [alldimIDs;fg.v[fid].attributes["dimIDs"]]
    end
    for cid in cond
      alldimIDs = [alldimIDs;fg.v[cid].attributes["data"].dimIDs]
      # alldimIDs = [alldimIDs;fg.v[cid].attributes["dimIDs"]]
    end

    for fid in frtl
        usefcts = []
        for fct in out_neighbors(fg.v[fid],fg.g)
            #println("Checking fct=$(fct.label)")
            if fct.attributes["data"].potentialused!=true #fct.attributes["potentialused"]!=true ## USED TO HAVE == AND continue end here
                if length(out_neighbors(fct, fg.g))==1
                    appendUseFcts!(usefcts, fg.IDs[out_neighbors(fct, fg.g)[1].label], fct, fid)
                    #usefcts = [usefcts;(fg.IDs[out_neighbors(fct, fg.g)[1].label], fct, fid)]
                    fct.attributes["data"].potentialused = true #fct.attributes["potentialused"] = true
                end
                for sepSearch in out_neighbors(fct, fg.g)
                    if (fg.IDs[sepSearch.label] == fid)
                        continue # skip the fid itself
                    end
                    sea = findmin(abs(allids-fg.IDs[sepSearch.label]))
                    if sea[1]==0.0
                        appendUseFcts!(usefcts, fg.IDs[sepSearch.label], fct, fid)
                        # usefcts = [usefcts;(fg.IDs[sepSearch.label], fct, fid)]
                        fct.attributes["data"].potentialused = true #fct.attributes["potentialused"] = true
                    end
                end
            end
        end
        cliq.attributes["potentials"]=[cliq.attributes["potentials"];usefcts]
    end
    return nothing
end

function getCliquePotentials!(fg::FactorGraph, bt::BayesTree, chkcliq::Int64)
    getCliquePotentials!(fg, bt.cliques[chkcliq])
end

function cliqPotentialIDs(cliq::Graphs.ExVertex)
  potIDs = Int[]
  for idfct in cliq.attributes["potentials"]
    push!(potIDs,idfct[1])
  end
  return potIDs
end

function collectSeparators(bt::BayesTree, cliq::Graphs.ExVertex)
  allseps = Int[]
  for child in out_neighbors(cliq, bt.bt)#tree
    # if length(child.attributes["conditIDs"]) > 0
      allseps = [allseps; child.attributes["conditIDs"]]
    # end
  end
  return allseps
end

function compCliqAssocMatrices!(bt::BayesTree, cliq::Graphs.ExVertex)
  frtl = cliq.attributes["frontalIDs"]
  cond = cliq.attributes["conditIDs"]
  inmsgIDs = collectSeparators(bt, cliq)
  potIDs = cliqPotentialIDs(cliq)
  # Construct associations matrix here
  # matrix has variables are columns, and messages/constraints as rows
  cols = [frtl;cond]
  cliq.attributes["inmsgIDs"] = inmsgIDs
  cliq.attributes["potIDs"] = potIDs
  cliqAssocMat = Array{Bool,2}(length(potIDs), length(cols))
  cliqMsgMat = Array{Bool,2}(length(inmsgIDs), length(cols))
  fill!(cliqAssocMat, false)
  fill!(cliqMsgMat, false)
  for j in 1:length(cols)
    for i in 1:length(inmsgIDs)
      if cols[j] == inmsgIDs[i]
        cliqMsgMat[i,j] = true
      end
    end
    for i in 1:length(potIDs)
      idfct = cliq.attributes["potentials"][i]
      if idfct[2].index == potIDs[i] # sanity check on clique potentials ordering
        for vertidx in idfct[2].attributes["data"].fncargvID #for vert in idfct[2].attributes["fnc"].Xi
          if vertidx == cols[j]
            cliqAssocMat[i,j] = true
          end
        end
        # for vert in idfct[2].attributes["data"].fnc.Xi #for vert in idfct[2].attributes["fnc"].Xi
        #   if vert.index == cols[j]
        #     cliqAssocMat[i,j] = true
        #   end
        # end
      else
        prtslperr("compCliqAssocMatrices! -- potential ID ordering was lost")
      end
    end
  end
  cliq.attributes["cliqAssocMat"] = cliqAssocMat
  cliq.attributes["cliqMsgMat"] = cliqMsgMat
  nothing
end

function getCliqMat(cliq::Graphs.ExVertex; showmsg=true)
  assocMat = cliq.attributes["cliqAssocMat"]
  msgMat = cliq.attributes["cliqMsgMat"]
  mat = showmsg ? [assocMat;msgMat] : assocMat
  # @show size(mat), size(assocMat), size(msgMat)
  return mat
end

function spyCliqMat(cliq::Graphs.ExVertex; showmsg=true)
  mat = getCliqMat(cliq, showmsg=showmsg)
  Gadfly.spy(mat)
end

function countSkips(bt::BayesTree)
  skps = 0
  for cliq in bt.cliques
    m = getCliqMat(cliq[2])
    mi = map(Int,m)
    skps += sum(map(Int,sum(mi,1) .== 1))
  end
  return skps
end

function skipThroughMsgsIDs(cliq::Graphs.ExVertex)
  numfrtl1 = floor(Int,length(cliq.attributes["frontalIDs"])+1)
  condAssocMat = cliq.attributes["cliqAssocMat"][:,numfrtl1:end]
  condMsgMat = cliq.attributes["cliqMsgMat"][:,numfrtl1:end]
  mat = [condAssocMat;condMsgMat];
  mab = sum(map(Int,mat),1) .== 1
  mabM = sum(map(Int,condMsgMat),1) .== 1
  mab = mab & mabM
  # rang = 1:size(condMsgMat,2)
  msgidx = cliq.attributes["conditIDs"][collect(mab)]
  return msgidx
end

function directAssignmentIDs(cliq::Graphs.ExVertex)
  assocMat = cliq.attributes["cliqAssocMat"]
  msgMat = cliq.attributes["cliqMsgMat"]
  mat = [assocMat;msgMat];
  mab = sum(map(Int,mat),1) .== 1
  mabA = sum(map(Int,assocMat),1) .== 1
  mab = mab & mabA
  frtl = cliq.attributes["frontalIDs"]
  cond = cliq.attributes["conditIDs"]
  cols = [frtl;cond]
  return cols[collect(mab)]
end

function directFrtlMsgIDs(cliq::Graphs.ExVertex)
  numfrtl = length(cliq.attributes["frontalIDs"])
  frntAssocMat = cliq.attributes["cliqAssocMat"][:,1:numfrtl]
  frtlMsgMat = cliq.attributes["cliqMsgMat"][:,1:numfrtl]
  mat = [frntAssocMat; frtlMsgMat];
  mab = sum(map(Int,mat),1) .== 1
  mabM = sum(map(Int,frtlMsgMat),1) .== 1
  mab = mab & mabM
  return cliq.attributes["frontalIDs"][collect(mab)]
end

function mcmcIterationIDs(cliq::Graphs.ExVertex)
  assocMat = cliq.attributes["cliqAssocMat"]
  msgMat = cliq.attributes["cliqMsgMat"]
  mat = [assocMat;msgMat];
  sum(sum(map(Int,mat),1)) == 0 ? error("mcmcIterationIDs -- unaccounted variables") : nothing
  mab = 1 .< sum(map(Int,mat),1)
  frtl = cliq.attributes["frontalIDs"]
  cond = cliq.attributes["conditIDs"]
  cols = [frtl;cond]
  return cols[collect(mab)]
end

function setCliqMCIDs!(cliq::Graphs.ExVertex)
  cliq.attributes["itervarIDs"] = mcmcIterationIDs(cliq)
  cliq.attributes["msgskipIDs"] = skipThroughMsgsIDs(cliq)
  cliq.attributes["directvarIDs"] = directAssignmentIDs(cliq)
  cliq.attributes["directFrtlMsgIDs"] = directFrtlMsgIDs(cliq)
  nothing
end


# function drawCliqAssocMat(fgl::FactorGrapha, cliq::Graphs.ExVertex)
#   cliqAssocMat = cliq.attributes["cliqAssocMat"]
#   ro, co = size(cliqAssocMat)
#   [@show fg.f[i].label for i in tree.cliques[1].attributes["potIDs"]];
#
#   showmat = Array{ASCIIString
#
# end

# post order tree traversal and build potential functions
function buildCliquePotentials(fg::FactorGraph, bt::BayesTree, cliq::Graphs.ExVertex)
    for child in out_neighbors(cliq, bt.bt)#tree
        buildCliquePotentials(fg, bt, child)
    end
    println("Get potentials $(cliq.attributes["label"])");
    getCliquePotentials!(fg, bt, cliq);

    compCliqAssocMatrices!(bt, cliq);
    setCliqMCIDs!(cliq);

    return Union{}
end
