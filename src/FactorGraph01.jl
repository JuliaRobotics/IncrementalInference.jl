
reshapeVec2Mat(vec::Vector, rows::Int) = reshape(vec, rows, round(Int,length(vec)/rows))
# function reshapeVec2Mat(vec::Vector, rows::Int)
#   M = reshape(vec, rows, round(Int,length(vec)/rows))
#   return ndims(M) < 2 ? (M')' : M
# end



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
  vert = dlapi.addvertex!(fg.g, ExVertex(fg.id,lbl))
  dlapi.setupvertgraph!(fg, vert, fg.id, lbl) #fg.v[fg.id]
  # fg.v[fg.id] = Graphs.add_vertex!(fg.g, ExVertex(fg.id,lbl))
  # fg.IDs[lbl] = fg.id # TODO -- API issue
  # fg.v[fg.id].attributes = Graphs.AttributeDict()
  # fg.v[fg.id].attributes["label"] = lbl
  dodims = fg.dimID+1

  # TODO -- vert should not lose information here
  setDefaultNodeData!(vert, initval, stdev, dodims, N) #fg.v[fg.id]

  fg.dimID+=size(initval,1) # rows indicate dimensions, move to last dimension
  push!(fg.nodeIDs,fg.id)
  return vert #fg.v[fg.id]
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
    vert = dlapi.getvertex(fg,p)
    if vert.attributes["data"].BayesNetVertID == 0   #fg.v[p].
      fg.bnid+=1
      # fg.bnverts[p] = Graphs.add_vertex!(fg.bn, ExVertex(fg.bnid,string("BayesNet",fg.bnid)))
      vert.attributes["data"].BayesNetVertID = p # fg.v[p] ##fg.bnverts[p]
      # fg.bnverts[p].attributes["label"] = fg.v[p].attributes["label"]
    else
      println("addBayesNetVerts -- something is very wrong, should not have a Bayes net vertex")
    end
  end
end

function addConditional!(fg::FactorGraph, vertID, lbl, Si)
  bnv = dlapi.getvertex(fg, vertID) #fg.v[vertID]
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
    push!(Xi, dlapi.getvertex(fg, s))  # push!(Xi,fg.v[s])
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
      for fct in out_neighbors(dlapi.getvertex(fg,v),fg.g) # (fg.v[v]
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
      dlapi.getvertex(fg,v).attributes["data"].eliminated = true # fg.v[v].
      # fg.v[v].attributes["eliminated"] = true

      #add marginal on remaining variables... ? f(xyz) = f(x | yz) f(yz)
      # new function between all Si
      addChainRuleMarginal!(fg, Si)
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
  return kde!(getVal(v), "lcv") # TODO - getBW(v)
end

function getVertKDE(fgl::FactorGraph, id::Int64)
  v = dlapi.getvertex(fgl,id)  # v = fgl.v[fgl.IDs[lbl]]
  return getKDE(v)
end
function getVertKDE(fgl::FactorGraph, lbl::ASCIIString)
  v = dlapi.getvertex(fgl,lbl)  # v = fgl.v[fgl.IDs[lbl]]
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
