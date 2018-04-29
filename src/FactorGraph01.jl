
reshapeVec2Mat(vec::Vector, rows::Int) = reshape(vec, rows, round(Int,length(vec)/rows))
# function reshapeVec2Mat(vec::Vector, rows::Int)
#   M = reshape(vec, rows, round(Int,length(vec)/rows))
#   return ndims(M) < 2 ? (M')' : M
# end

# get vertex from factor graph according to label symbol "x1"
getVert(fgl::FactorGraph, lbl::Symbol; api::DataLayerAPI=dlapi, nt::Symbol=:var) = api.getvertex(fgl, lbl, nt=nt)
getVert(fgl::FactorGraph, id::Int; api::DataLayerAPI=dlapi) = api.getvertex(fgl, id)


# TODO -- upgrade to dedicated memory location in Graphs.jl
# see JuliaArchive/Graphs.jl#233
getData(v::Graphs.ExVertex) = v.attributes["data"]
# Convenience functions
getData(fgl::FactorGraph, lbl::Symbol; api::DataLayerAPI=dlapi) = getData(getVert(fgl, lbl, api=api))
getData(fgl::FactorGraph, id::Int; api::DataLayerAPI=dlapi) = getData(getVert(fgl, id, api=api))

function setData!(v::Graphs.ExVertex, data)
  v.attributes["data"] = data
  nothing
end

function getVal(v::Graphs.ExVertex)
  return getData(v).val
end
function getVal(v::Graphs.ExVertex, idx::Int)
  return getData(v).val[:,idx]
end

# Convenience function to get values for given variable label
function getVal(fgl::FactorGraph, lbl::Symbol; api::DataLayerAPI=dlapi)
  #getVal(dlapi.getvertex(fgl, lbl))
  getVal(getVert(fgl, lbl, api=api))
end
function getVal(fgl::FactorGraph, exvertid::Int; api::DataLayerAPI=dlapi)
  # getVal(dlapi.getvertex(fgl, exvertid))
  getVal(getVert(fgl, exvertid, api=api))
end

function getfnctype(vertl::Graphs.ExVertex)
  data = getData(vertl)
  if typeof(data).name.name == :VariableNodeData
    return VariableNodeData
  end
  data.fnc.usrfnc!
end

function getfnctype(fgl::FactorGraph, exvertid::Int; api::DataLayerAPI=dlapi)
  #
  # data = getData(fgl, exvertid, api=api)
  # data.fnc.usrfnc!
  getfnctype(getVert(fgl, exvertid, api=api))
end

# setVal! assumes you will update values to database separate, this used for local graph mods only
function setVal!(v::Graphs.ExVertex, val::Array{Float64,2})
  v.attributes["data"].val = val
  nothing
end
function getBWVal(v::Graphs.ExVertex)
  return getData(v).bw
end
function setBW!(v::Graphs.ExVertex, bw::Array{Float64,2})
  v.attributes["data"].bw = bw
  nothing
end
function setVal!(v::Graphs.ExVertex, val::Array{Float64,2}, bw::Array{Float64,2})
  setVal!(v,val)
  setBW!(v,bw)
  nothing
end
function setVal!(v::Graphs.ExVertex, val::Array{Float64,2}, bw::Vector{Float64})
  setVal!(v,val,reshape(bw,length(bw),1)) #(bw')')
  nothing
end
function setValKDE!(v::Graphs.ExVertex, val::Array{Float64,2})
  p = kde!(val)
  setVal!(v,val,getBW(p)[:,1]) # TODO -- this can be little faster
  nothing
end
function setValKDE!(v::Graphs.ExVertex, em::EasyMessage)
  setVal!(v, em.pts, em.bws ) # getBW(p)[:,1]
  nothing
end
function setValKDE!(v::Graphs.ExVertex, p::BallTreeDensity)
  pts = getPoints(p)
  setVal!(v, pts, getBW(p)[:,1]) # BUG ...al!(., val, . ) ## TODO -- this can be little faster
  nothing
end
setVal!(v::Graphs.ExVertex, em::EasyMessage) = setValKDE!(v, em)
setVal!(v::Graphs.ExVertex, p::BallTreeDensity) = setValKDE!(v, p)


function kde!(em::EasyMessage)
  return kde!(em.pts,em.bws)
end


# TODO -- there should be a better way, without retrieving full vertex
getOutNeighbors(fgl::FactorGraph, v::ExVertex; api::DataLayerAPI=dlapi, needdata::Bool=false, ready::Int=1,backendset::Int=1 ) = api.outneighbors(fgl, v, needdata=needdata, ready=ready, backendset=backendset )
getOutNeighbors(fgl::FactorGraph, vertid::Int; api::DataLayerAPI=dlapi, needdata::Bool=false, ready::Int=1,backendset::Int=1 ) = api.outneighbors(fgl, api.getvertex(fgl,vertid), needdata=needdata, ready=ready, backendset=backendset )

function updateFullVert!(fgl::FactorGraph, exvert::ExVertex;
            api::DataLayerAPI=IncrementalInference.dlapi,
            updateMAPest::Bool=false  )
  #
  warn("use of updateFullVert! should be clarified for local or remote operations.")
  api.updatevertex!(fgl, exvert, updateMAPest=updateMAPest)
end


function setDefaultNodeData!(v::Graphs.ExVertex, initval::Array{Float64,2},
                             stdev::Array{Float64,2}, dodims::Int, N::Int, dims::Int;
                             gt=nothing, initialized::Bool=true,
                             softtype=nothing)
  # TODO review and refactor this function, exists as legacy from pre-v0.3.0
  # this should be the only function allocating memory for the node points (unless number of points are changed)
  data = nothing
  if initialized
    if size(initval,2) < N && size(initval, 1) == dims
      warn("setDefaultNodeData! -- deprecated use of stdev.")
      p = kde!(initval,diag(stdev));
      pN = resample(p,N)
    elseif size(initval,2) < N && size(initval, 1) != dims
      println("Node value memory allocated but not initialized")
      pN = kde!(randn(dims, N));
    else
      pN = kde!(initval)
    end
    # dims = size(initval,1) # rows indicate dimensions
    sp = round.(Int,linspace(dodims,dodims+dims-1,dims))
    gbw = getBW(pN)[:,1]
    gbw2 = Array{Float64}(length(gbw),1)
    gbw2[:,1] = gbw[:]
    pNpts = getPoints(pN)
    data = VariableNodeData(initval, stdev, pNpts,
                            gbw2, Int[], sp,
                            dims, false, 0, Int[], gt, softtype, true) #initialized
  else
      sp = round.(Int,linspace(dodims,dodims+dims-1,dims))
      data = VariableNodeData(initval, stdev, zeros(dims, N),
                              zeros(dims,1), Int[], sp,
                              dims, false, 0, Int[], gt, softtype, false) #initialized
  end
  #
  v.attributes["data"] = data

  nothing
end

function addNewVarVertInGraph!(fgl::FactorGraph, vert::Graphs.ExVertex, id::Int, lbl::Symbol, ready::Int)
  vert.attributes = Graphs.AttributeDict() #fg.v[fg.id]
  vert.attributes["label"] = string(lbl) #fg.v[fg.id]
  fgl.IDs[lbl] = id

  # used for cloudgraph solving
  vert.attributes["ready"] = ready
  vert.attributes["backendset"] = 0

  # fgl.g.vertices[id] = vert # will be inserted during addvertex! call
  # fgl.v[id] = vert #  -- this is likely not required, but is used in subgraph methods
  nothing
end

# must set either dims or initval for proper initialization
# Add node to graph, given graph struct, labal, init values,
# std dev [TODO -- generalize], particle size and ready flag for concurrency
function addNode!(fg::FactorGraph,
      lbl::Symbol,
      initval::Array{Float64}=zeros(1,1),
      stdev::Array{Float64}=ones(1,1); # this is bad and should be removed TODO
      N::Int=100,
      ready::Int=1,
      labels::Vector{T}=String[],
      api::DataLayerAPI=dlapi,
      uid::Int=-1,
      dims::Int=-1  ) where {T <: AbstractString}
  #
  warn("this addNode! will be deprecated, please use FactorGraph01.jl:addNode!(fg::FactorGraph, lbl::Symbol, softtype::Type{T}).")
  currid = fg.id+1
  if uid==-1
    fg.id=currid
  else
    currid = uid
  end
  dims = dims != -1 ? dims : size(initval,1)

  lblstr = string(lbl)
  vert = ExVertex(currid,lblstr)
  addNewVarVertInGraph!(fg, vert, currid, lbl, ready)
  # dlapi.setupvertgraph!(fg, vert, currid, lbl) #fg.v[currid]
  dodims = fg.dimID+1
  # TODO -- vert should not loose information here
  setDefaultNodeData!(vert, initval, stdev, dodims, N, dims) #fg.v[currid]

  vnlbls = deepcopy(labels)
  push!(vnlbls, fg.sessionname)
  # addvert!(fg, vert, api=api)
  api.addvertex!(fg, vert, labels=vnlbls) #fg.g ##vertr =

  fg.dimID+=dims # rows indicate dimensions, move to last dimension
  push!(fg.nodeIDs, currid)

  return vert #fg.v[fg.id]
end

"""
$(SIGNATURES)

Add a node (variable) to a graph. Use this over the other dispatches.
"""
function addNode!(fg::FactorGraph,
      lbl::Symbol,
      softtype::Type{<:InferenceVariable};
      N::Int=100,
      autoinit=true,  # does init need to be separate from ready? TODO
      ready::Int=1,
      labels::Vector{<:AbstractString}=String[],
      api::DataLayerAPI=dlapi,
      uid::Int=-1,
      dims::Int=-1  ) # where {T , S }
  #
  currid = fg.id+1
  if uid==-1
    fg.id=currid
  else
    currid = uid
  end
  st = softtype()
  dims = dims != -1 ? dims : st.dims

  lblstr = string(lbl)
  vert = ExVertex(currid,lblstr)
  addNewVarVertInGraph!(fg, vert, currid, lbl, ready)
  # dlapi.setupvertgraph!(fg, vert, currid, lbl) #fg.v[currid]
  dodims = fg.dimID+1
  setDefaultNodeData!(vert, zeros(dims,N), zeros(0,0), dodims, N, dims, initialized=!autoinit, softtype=st) #fg.v[currid]

  vnlbls = union(string.(labels), st.labels)
  push!(vnlbls, fg.sessionname)
  # addvert!(fg, vert, api=api)
  api.addvertex!(fg, vert, labels=vnlbls) #fg.g ##vertr =

  fg.dimID+=dims # rows indicate dimensions, move to last dimension
  push!(fg.nodeIDs, currid)

  vert
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
  cols = Array{Int,1}()
  push!(cols,0)
  rows = Array{Int,1}()
  for v in vA
      push!(vals, getVal(v))
      c = size(vals[end],2)
      r = size(vals[end],1)
      push!(cols, floor(Int,c))
      push!(rows, floor(Int,r))
  end
  cols = cumsum(cols)
  sc = cols[end]
  rw = floor(Int,rows[1])
  val = Array{Float64,2}(rw, sc)
  for i in 1:(len-1)
      val[:,(cols[i]+1):cols[i+1]] = vals[i]
  end
  val[:,(cols[len]+1):cols[len+1]] = vals[len] # and the last one
  return val
end

function registerCallback!(fgl::FactorGraph, fnc::Function)
  m = Symbol(typeof(fnc).name.module)
  fgl.registeredModuleFunctions[m] = fnc
  nothing
end

# idea for registering, will likely be removed again
# function createregistercallback!(fnc::Function)
#   fgl = initfg()
#   m = Symbol(typeof(fnc).name.module)
#   fgl.registeredModuleFunctions[m] = fnc
#   nothing
# end
#
function prepareparamsarray!(ARR::Array{Array{Float64,2},1},Xi::Vector{Graphs.ExVertex}, N::Int, solvefor::Int)
  LEN = Int[]
  maxlen = N
  count = 0
  sfidx = 0
  for xi in Xi
    push!(ARR, getVal(xi))
    len = size(ARR[end], 2)
    push!(LEN, len)
    if len > maxlen
      maxlen = len
    end
    count += 1
    if xi.index == solvefor
      sfidx = count #xi.index
    end
  end
  SAMP=LEN.<maxlen
  for i in 1:count
    if SAMP[i]
      ARR[i] = KernelDensityEstimate.sample(getKDE(Xi[i]), maxlen)[1]
    end
  end
  # we are generating a proposal distribution, not direct replacement for existing
  if sfidx > 0 ARR[sfidx] = deepcopy(ARR[sfidx]) end
  return maxlen, sfidx
end

function prepgenericwrapper{T <: FunctorInferenceType}(
      Xi::Vector{Graphs.ExVertex},
      usrfnc::T,
      samplefnc::Function )
  #
  ARR = Array{Array{Float64,2},1}()
  maxlen, sfidx = prepareparamsarray!(ARR, Xi, 0, 0)
  # test if specific zDim or partial constraint used
  fldnms = fieldnames(usrfnc)
  # sum(fldnms .== :zDim) >= 1
  return GenericWrapParam{T}(usrfnc, ARR, 1, 1, (zeros(0,1),), samplefnc, sum(fldnms .== :zDim) >= 1, sum(fldnms .== :partial) >= 1)
end

function setDefaultFactorNode!{T <: Union{FunctorInferenceType, InferenceType}}(
      fgl::FactorGraph,
      vert::Graphs.ExVertex,
      Xi::Vector{Graphs.ExVertex},
      usrfnc::T  )
  #
  ftyp = typeof(usrfnc) # maybe this can be T
  m = Symbol(ftyp.name.module)
  # samplefnc2 = fgl.registeredModuleFunctions[m]
  gwpf = prepgenericwrapper(Xi, usrfnc, getSample) #samplefnc2)

  data = FunctionNodeData{GenericWrapParam{T}}(Int[], false, false, Int[], m, gwpf)
  vert.attributes["data"] = data

  nothing
end

function addNewFncVertInGraph!(fgl::FactorGraph, vert::Graphs.ExVertex, id::Int, lbl::Symbol, ready::Int)
  vert.attributes = Graphs.AttributeDict() #fg.v[fg.id]
  vert.attributes["label"] = lbl #fg.v[fg.id]
  # fgl.f[id] = vert #  -- not sure if this is required, using fg.g.vertices
  fgl.fIDs[lbl] = id # fg.id

  # used for cloudgraph solving
  vert.attributes["ready"] = ready
  vert.attributes["backendset"] = 0

  # for graphviz drawing
  vert.attributes["shape"] = "point"
  vert.attributes["width"] = 0.2
  nothing
end
addNewFncVertInGraph!{T <: AbstractString}(fgl::FactorGraph, vert::Graphs.ExVertex, id::Int, lbl::T, ready::Int) =
    addNewFncVertInGraph!(fgl,vert, id, Symbol(lbl), ready)

function isInitialized(vert::Graphs.ExVertex)::Bool
  return getData(vert).initialized
end
function isInitialized(fgl::FactorGraph, vsym::Symbol)::Bool
  # TODO, make cloudgraphs work and make faster by avoiding all the getVerts
  isInitialized(getVert(fgl, vsym))
end

"""
    doautoinit!(fg, Xi[,api=dlapi])

initialize destination variable nodes based on this factor in factor graph, fg, generally called
during addFactor!. Destination factor is first (singletons) or second (dim 2 pairwise) variable vertex in Xi.
"""
function doautoinit!(fgl::FactorGraph, Xi::Vector{Graphs.ExVertex}; api::DataLayerAPI=dlapi, singles::Bool=true)
  # Mighty inefficient function, since we only need very select fields nearby from a few neighboring nodes
  # do double depth search for variable nodes
  # TODO this should maybe stay localapi only...
  for xi in Xi
    if !isInitialized(xi)
       # @show "doautoinit!", size(getVal(xi))
      vsym = Symbol(xi.label)
      neinodes = ls(fgl, vsym)
      if (length(neinodes) > 1 || singles) # && !isInitialized(xi)
        # Which of the factors can be used for initialization
        useinitfct = Symbol[]
        # println("Consider all pairwise factors connected to $vsym...")
        for xifct in neinodes #potntlfcts
          xfneivarnodes = lsf(fgl, xifct)
          for vsym2 in xfneivarnodes
            # println("find all variables that are initialized for $vsym2")
            vert2 = getVert(fgl, vsym2)
            if (isInitialized(vert2) && sum(useinitfct .== xifct) == 0 ) || length(xfneivarnodes) == 1      # OR singleton  TODO get faster version of isInitialized for database version
              # println("adding $xifct to init factors list")
              push!(useinitfct, xifct)
            end
          end
        end
        # println("Consider all singleton (unary) factors to $vsym...")

        # calculate the predicted belief over $vsym
        pts = predictbelief(fgl, vsym, useinitfct, api=api)
        # println("doautoinit! Past predictbelief...")
        setValKDE!(xi, pts)
        getData(xi).initialized = true
        # println("doautoinit! just before update vertex...")
        api.updatevertex!(fgl, xi, updateMAPest=false)
      end
    end
  end

  # len = length(Xi)
  # # pts = Array{Float64,2}()
  # if length(Xi) == 1
  #   pts = evalFactor2(fgl, fc, Xi[1].index)
  #   setValKDE!(Xi[1], pts)
  #   api.updatevertex!(fgl, Xi[1], updateMAPest=false)
  # elseif len == 2
  #   pts = evalFactor2(fgl, fc, Xi[2].index)
  #   setValKDE!(Xi[2], pts)
  #   api.updatevertex!(fgl, Xi[2], updateMAPest=false)
  # else
  #   # consider specifying an init order in the constraint type
  #   # also consider taking product between all incoming densities which have been inited
  #   error("don't know how to autoinit with pairwise dimension > 2")
  # end
  # println("doautoinit! Done in autoinit!")
  nothing
end

function ensureAllInitialized!(fgl::FactorGraph; api::DataLayerAPI=dlapi)
  xx, xl = ls(fgl)
  allvarnodes = union(xx, xl)
  for vsym in allvarnodes
    if !isInitialized(fgl, vsym)
      println("$vsym is not initialized, and will do so now...")
      doautoinit!(fgl, Graphs.ExVertex[getVert(fgl, vsym, api=api);], api=api, singles=true)
    end
  end
  nothing
end

function assembleFactorName(fgl::FactorGraph, Xi::Vector{Graphs.ExVertex}; maxparallel::Int=50)
  namestring = ""
  for vert in Xi #f.Xi
    namestring = string(namestring,vert.attributes["label"])
  end
  for i in 1:maxparallel
    tempnm = string(namestring, "f$i")
    if !haskey(fgl.fIDs, Symbol(tempnm))
      namestring = tempnm
      break
    end
    i != maxparallel ? nothing : error("Cannot currently add more than $(maxparallel) factors in parallel, please open an issue if this is too restrictive.")
  end
  return namestring
end

function addFactor!(fgl::FactorGraph,
      Xi::Vector{Graphs.ExVertex},
      usrfnc::R;
      ready::Int=1,
      api::DataLayerAPI=dlapi,
      labels::Vector{T}=String[],
      uid::Int=-1,
      autoinit::Bool=true) where
        {R <: Union{FunctorInferenceType, InferenceType},
         T <: AbstractString}
  #
  currid = fgl.id+1
  if uid==-1
    fgl.id=currid
  else
    currid = uid
  end
  namestring = assembleFactorName(fgl, Xi)
  # fgl.id+=1
  newvert = ExVertex(currid,namestring)
  addNewFncVertInGraph!(fgl, newvert, currid, namestring, ready)
  setDefaultFactorNode!(fgl, newvert, Xi, deepcopy(usrfnc))
  push!(fgl.factorIDs,currid)

  for vert in Xi
    push!(newvert.attributes["data"].fncargvID, vert.index)
  end

  fnlbls = deepcopy(labels)
  push!(fnlbls, "FACTOR")
  push!(fnlbls, fgl.sessionname)
  # TODO -- multiple accesses to DB with this method, must refactor!
  newvert = api.addvertex!(fgl, newvert, labels=fnlbls)  # used to be two be three lines up ##fgl.g
  for vert in Xi
    api.makeaddedge!(fgl, vert, newvert)
  end

  # TODO change this operation to update a conditioning variable
  if autoinit
    doautoinit!(fgl, Xi, api=api, singles=false)
  end

  return newvert
end

function addFactor!(
      fgl::FactorGraph,
      xisyms::Vector{Symbol},
      usrfnc::R;
      ready::Int=1,
      api::DataLayerAPI=dlapi,
      labels::Vector{T}=String[],
      uid::Int=-1,
      autoinit::Bool=true ) where
        {R <: Union{FunctorInferenceType, InferenceType},
         T <: AbstractString}
  #
  verts = Vector{Graphs.ExVertex}()
  for xi in xisyms
      push!( verts, api.getvertex(fgl,xi) )
  end
  addFactor!(fgl, verts, usrfnc, ready=ready, api=api, labels=labels, uid=uid, autoinit=autoinit)
end


function prtslperr(s)
  println(s)
  sleep(0.1)
  error(s)
end

# for computing the Bayes Net-----------------------------------------------------
function getEliminationOrder(fg::FactorGraph; ordering::Symbol=:qr)
    s = fg.nodeIDs
    lens = length(s)
    sf = fg.factorIDs
    lensf = length(sf)
    adjm, dictpermu = adjacency_matrix(fg.g,returnpermutation=true)
    permuteds = Vector{Int}(lens)
    permutedsf = Vector{Int}(lensf)
    for j in 1:length(dictpermu)
      semap = 0
      for i in 1:lens
        if dictpermu[j] == s[i]
          permuteds[i] = j#dictpermu[j]
          semap += 1
          if semap >= 2  break; end
        end
      end
      for i in 1:lensf
        if dictpermu[j] == sf[i]
          permutedsf[i] = j#dictpermu[j]
          semap += 1
          if semap >= 2  break; end
        end
      end
    end

      A=convert(Array{Int},adjm[permutedsf,permuteds]) # TODO -- order seems brittle
      p = Int[]
      if ordering==:chol
        p = cholfact(A'A,:U,Val{true})[:p] #,pivot=true
      elseif ordering==:qr
        q,r,p = qr(A,Val{true})
      else
        prtslperr("getEliminationOrder -- cannot do the requested ordering $(ordering)")
      end

    # we need the IDs associated with the Graphs.jl and our Type fg
    return  dictpermu[permuteds[p]] # fg.nodeIDs[p]
end


# lets create all the vertices first and then deal with the elimination variables thereafter
function addBayesNetVerts!(fg::FactorGraph, elimOrder::Array{Int,1})
  for p in elimOrder
    vert = getVert(fg, p, api=localapi)
    @show vert.label, getData(vert).BayesNetVertID
    if getData(vert).BayesNetVertID == 0
      fg.bnid+=1
      vert.attributes["data"].BayesNetVertID = p
      localapi.updatevertex!(fg, vert)
    else
      println("addBayesNetVerts -- something is very wrong, should not have a Bayes net vertex")
    end
  end
end

function addConditional!(fg::FactorGraph, vertID::Int, lbl, Si)
  bnv = getVert(fg, vertID, api=localapi) #fg.v[vertID]
  bnvd = getData(bnv) # bnv.attributes["data"]
  bnvd.separator = Si
  for s in Si
    push!(bnvd.BayesNetOutVertIDs, s)
  end
  localapi.updatevertex!(fg, bnv)
  nothing
end

function addChainRuleMarginal!(fg::FactorGraph, Si)
  lbls = String[]
  genmarg = GenericMarginal()
  Xi = Graphs.ExVertex[]
  for s in Si
    push!(Xi, getVert(fg, s, api=localapi))
  end
  println("adding marginal to")
  for x in Xi @show x.index end
  addFactor!(fg, Xi, genmarg, api=localapi)
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
        i = 0
        for id in alleids
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

      # all factors adjacent to this variable
      fi = Int[]
      Si = Int[]
      gm = ExVertex[]
      # TODO -- optimize outneighbor calls like this
      vert = localapi.getvertex(fg,v)
      for fct in localapi.outneighbors(fg, vert)
        if (getData(fct).eliminated != true)
          push!(fi, fct.index)
          for sepNode in localapi.outneighbors(fg, fct)
            if sepNode.index != v && length(findin(sepNode.index,Si)) == 0
              push!(Si,sepNode.index)
            end
          end
          getData(fct).eliminated = true #fct.attributes["data"].eliminated = true
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

      # TODO -- remove links from current vertex to any marginals
      rmVarFromMarg(fg, vert, gm)

      #add marginal on remaining variables... ? f(xyz) = f(x | yz) f(yz)
      # new function between all Si
      addChainRuleMarginal!(fg, Si)

    end
    nothing
end

# some plotting functions on the factor graph
function stackVertXY(fg::FactorGraph, lbl::String)
    v = dlapi.getvertex(fg,lbl)
    vals = getVal(v)
    X=vec(vals[1,:])
    Y=vec(vals[2,:])
    return X,Y
end

function getKDE(v::Graphs.ExVertex)
  return kde!(getVal(v), getBWVal(v)[:,1])
end

function getVertKDE(v::Graphs.ExVertex)
  return getKDE(v)
end

function getVertKDE(fgl::FactorGraph, id::Int; api::DataLayerAPI=dlapi)
  v = api.getvertex(fgl,id)
  return getKDE(v)
end
function getVertKDE(fgl::FactorGraph, lbl::Symbol; api::DataLayerAPI=dlapi)
  v = api.getvertex(fgl,lbl)
  return getKDE(v)
end


function drawCopyFG(fgl::FactorGraph)
  fgd = deepcopy(fgl)
  for (sym,vid) in fgd.IDs
    delete!(fgd.g.vertices[vid].attributes,"data")
    !haskey(fgd.g.vertices[vid].attributes,"frtend") ? nothing : delete!(fgd.g.vertices[vid].attributes,"frtend")
  end
  for (sym,vid) in fgd.fIDs
    delete!(fgd.g.vertices[vid].attributes,"data")
    !haskey(fgd.g.vertices[vid].attributes,"frtend") ? nothing : delete!(fgd.g.vertices[vid].attributes,"frtend")
  end
  return fgd
end

# function writeGraphPdf(fgl::FactorGraph)
#   fgd = drawCopyFG(fgl)
#   println("Writing factor graph file")
#   fid = open("fg.dot","w+")
#   write(fid,Graphs.to_dot(fgd.g))
#   close(fid)
#   run(`dot fg.dot -Tpdf -o fg.pdf`)
#   nothing
# end




function expandEdgeListNeigh!(fgl::FactorGraph,
                              vertdict::Dict{Int,Graphs.ExVertex},
                              edgedict::Dict{Int,Graphs.Edge{Graphs.ExVertex}})
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
  edgedict::Dict{Int,Graphs.Edge{Graphs.ExVertex}},
  vertdict::Dict{Int,Graphs.ExVertex})

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
  edgedict = Dict{Int,Graphs.Edge{Graphs.ExVertex}}()
  for edge in edgelist
    edgedict[edge.index] = edge
  end
  return edgedict
end

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

# NOTICE, nodeIDs and factorIDs are not brough over by this method yet
# must sort out for incremental updates
function genSubgraph(fgl::FactorGraph,
    vertdict::Dict{Int,Graphs.ExVertex})
    # edgedict::Dict{Int,Graphs.Edge{Graphs.ExVertex}},

  fgseg = FactorGraph() # new handle for just a segment of the graph
  fgseg.g = Graphs.inclist(Graphs.ExVertex,is_directed=false)

  fgseg.v = Dict{Int,Graphs.ExVertex}()
  fgseg.f = Dict{Int,Graphs.ExVertex}()
  fgseg.IDs = Dict{AbstractString,Int}()
  fgseg.fIDs = Dict{AbstractString,Int}()

  # TODO -- convert to use empty constructor since Graphs.incdict now works
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

  vertdict = getShortestPathNeighbors(fgl, from=from, to=to, neighbors=neighbors)
  return genSubgraph(fgl, vertdict)
end

# explore all shortest paths combinations in verts, add neighbors and reference subgraph
function subgraphFromVerts(fgl::FactorGraph,
    verts::Dict{Int,Graphs.ExVertex};
    neighbors::Int=0  )

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

# explore all shortest paths combinations in verts, add neighbors and reference subgraph
# Using unique index into graph data structure
function subgraphFromVerts(fgl::FactorGraph,
    verts::Array{String,1};
    neighbors::Int=0  )

  vertdict = Dict{Int,Graphs.ExVertex}()
  for vert in verts
    id = -1
    if haskey(fgl.IDs, vert)
      id = fgl.IDs[Symbol(vert)]
    else
      error("FactorGraph01 only supports variable node subgraph search at this point")
    end
    vertdict[id] = fgl.g.vertices[id]
  end

  return subgraphFromVerts(fgl,vertdict,neighbors=neighbors)
end
