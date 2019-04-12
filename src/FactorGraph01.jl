
reshapeVec2Mat(vec::Vector, rows::Int) = reshape(vec, rows, round(Int,length(vec)/rows))
# function reshapeVec2Mat(vec::Vector, rows::Int)
#   M = reshape(vec, rows, round(Int,length(vec)/rows))
#   return ndims(M) < 2 ? (M')' : M
# end

# get vertex from factor graph according to label symbol "x1"
getVert(fgl::FactorGraph, lbl::Symbol; api::DataLayerAPI=dlapi, nt::Symbol=:var) = api.getvertex(fgl, lbl, nt=nt)
getVert(fgl::FactorGraph, id::Int; api::DataLayerAPI=dlapi) = api.getvertex(fgl, id)

"""
    $(TYPEDSIGNATURES)

Get the `::Symbol` name for a node with `id::Int`.
"""
getSym(fgl::FactorGraph, id::Int) = Symbol(fgl.g.vertices[id].label)

"""
    $SIGNATURES

Retrieve data structure stored in a node.
"""
getData(v::Graphs.ExVertex) = v.attributes["data"]
getData(fgl::FactorGraph, lbl::Symbol; api::DataLayerAPI=dlapi, nt=:var) = getData(getVert(fgl, lbl, api=api, nt=nt))
getData(fgl::FactorGraph, id::Int; api::DataLayerAPI=dlapi) = getData(getVert(fgl, id, api=api))
# TODO -- upgrade to dedicated memory location in Graphs.jl
# see JuliaArchive/Graphs.jl#233

function setData!(v::Graphs.ExVertex, data)
  v.attributes["data"] = data
  nothing
end

"""
   $(SIGNATURES)

Variable nodes softtype information holding a variety of meta data associated with the type of variable stored in that node of the factor graph.
"""
function getSofttype(vnd::VariableNodeData)
  return vnd.softtype
end
function getSofttype(v::ExVertex)
  getSofttype(getData(v))
end

"""
    $(SIGNATURES)

Convenience function to get point values sampled i.i.d from marginal of `lbl` variable in the current factor graph.
"""
getVal(vnd::VariableNodeData) = vnd.val
getVal(vnd::VariableNodeData, idx::Int) = vnd.val[:,idx]
function getVal(v::Graphs.ExVertex)
  return getVal(getData(v))
end
function getVal(v::Graphs.ExVertex, idx::Int)
  return getVal(getData(v),idx)
end
function getVal(fgl::FactorGraph, lbl::Symbol; api::DataLayerAPI=dlapi)
  #getVal(dlapi.getvertex(fgl, lbl))
  getVal(getVert(fgl, lbl, api=api))
end
function getVal(fgl::FactorGraph, exvertid::Int; api::DataLayerAPI=dlapi)
  # getVal(dlapi.getvertex(fgl, exvertid))
  getVal(getVert(fgl, exvertid, api=api))
end

"""
    $(SIGNATURES)

Get the number of points used for the current marginal belief estimate represtation for a particular variable in the factor graph.
"""
function getNumPts(v::Graphs.ExVertex)
  return size(getData(v).val,2)
end

function getfnctype(data)
  if typeof(data).name.name == :VariableNodeData
    return VariableNodeData
  end
  data.fnc.usrfnc!
end

function getfnctype(vertl::Graphs.ExVertex)
  data = getData(vertl)
  getfnctype(data)
end

function getfnctype(fgl::FactorGraph, exvertid::Int; api::DataLayerAPI=dlapi)
  #
  # data = getData(fgl, exvertid, api=api)
  # data.fnc.usrfnc!
  getfnctype(getVert(fgl, exvertid, api=api))
end

function getBW(vnd::VariableNodeData)
  return vnd.bw
end

# setVal! assumes you will update values to database separate, this used for local graph mods only
function setVal!(v::Graphs.ExVertex, val::Array{Float64,2})
  getData(v).val = val
  # v.attributes["data"].val = val
  nothing
end
function getBWVal(v::Graphs.ExVertex)
  return getData(v).bw
end
function setBW!(v::Graphs.ExVertex, bw::Array{Float64,2})
  getData(v).bw = bw
  # v.attributes["data"].bw = bw
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
function setVal!(fg::FactorGraph, sym::Symbol, val::Array{Float64,2}; api::DataLayerAPI=localapi)
  setVal!(api.getvertex(fg, sym), val)
end

"""
    $SIGNATURES

Set the point centers and bandwidth parameters of a variable node, also set `isInitialized=true` if `setinit::Bool=true` (as per default).
"""
function setValKDE!(v::Graphs.ExVertex, val::Array{Float64,2}, setinit::Bool=true)
  # recover softtype information
  sty = getSofttype(v)
  # @show sty.manifolds
  #
  p = AMP.manikde!(val, sty.manifolds)
  setVal!(v,val,getBW(p)[:,1]) # TODO -- this can be little faster
  setinit ? (getData(v).initialized = true) : nothing
  nothing
end
function setValKDE!(v::Graphs.ExVertex, em::EasyMessage, setinit::Bool=true)
  setVal!(v, em.pts, em.bws ) # getBW(p)[:,1]
  setinit ? (getData(v).initialized = true) : nothing
  nothing
end
function setValKDE!(v::Graphs.ExVertex, p::BallTreeDensity, setinit::Bool=true)
  pts = getPoints(p)
  setVal!(v, pts, getBW(p)[:,1]) # BUG ...al!(., val, . ) ## TODO -- this can be little faster
  setinit ? (getData(v).initialized = true) : nothing
  nothing
end
function setValKDE!(fgl::FactorGraph, sym::Symbol, p::BallTreeDensity; api::DataLayerAPI=dlapi, setinit::Bool=true)
  setValKDE!(getVert(fgl, sym, api=api), p, setinit)
  nothing
end
setVal!(v::Graphs.ExVertex, em::EasyMessage) = setValKDE!(v, em)
setVal!(v::Graphs.ExVertex, p::BallTreeDensity) = setValKDE!(v, p)

"""
    $(SIGNATURES)

Construct a BallTreeDensity KDE object from an IIF.EasyMessage object.
"""
function kde!(em::EasyMessage)
  return AMP.manikde!(em.pts,em.bws, em.manifolds)
end


# TODO -- there should be a better way, without retrieving full vertex
getOutNeighbors(fgl::FactorGraph, v::ExVertex; api::DataLayerAPI=dlapi, needdata::Bool=false, ready::Int=1,backendset::Int=1 ) = api.outneighbors(fgl, v, needdata=needdata, ready=ready, backendset=backendset )
getOutNeighbors(fgl::FactorGraph, vertid::Int; api::DataLayerAPI=dlapi, needdata::Bool=false, ready::Int=1,backendset::Int=1 ) = api.outneighbors(fgl, api.getvertex(fgl,vertid), needdata=needdata, ready=ready, backendset=backendset )

function updateFullVert!(fgl::FactorGraph, exvert::ExVertex;
            api::DataLayerAPI=IncrementalInference.dlapi,
            updateMAPest::Bool=false  )
  #
  @warn "use of updateFullVert! should be clarified for local or remote operations."
  api.updatevertex!(fgl, exvert, updateMAPest=updateMAPest)
end


function setDefaultNodeData!(v::Graphs.ExVertex,
                             dodims::Int,
                             N::Int,
                             dims::Int;
                             gt=nothing,
                             initialized::Bool=true,
                             dontmargin::Bool=false,
                             softtype=nothing)
  # TODO review and refactor this function, exists as legacy from pre-v0.3.0
  # this should be the only function allocating memory for the node points (unless number of points are changed)
  data = nothing
  if initialized
    # if size(initval,2) < N && size(initval, 1) == dims
    #   @warn "setDefaultNodeData! -- deprecated use of stdev."
    #   p = AMP.manikde!(initval,diag(stdev), softtype.manifolds);
    #   pN = resample(p,N)
    # if size(initval,2) < N && size(initval, 1) != dims
      @info "Node value memory allocated but not initialized"
      pN = AMP.manikde!(randn(dims, N), softtype.manifolds);
    # else
    #   pN = AMP.manikde!(initval, softtype.manifolds)
    # end
    # dims = size(initval,1) # rows indicate dimensions
    sp = Int[0;] #round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
    gbw = getBW(pN)[:,1]
    gbw2 = Array{Float64}(undef, length(gbw),1)
    gbw2[:,1] = gbw[:]
    pNpts = getPoints(pN)
    #initval, stdev
    data = VariableNodeData(pNpts,
                            gbw2, Int[], sp,
                            dims, false, 0, Int[], gt, softtype, true, false, dontmargin) #initialized
  else
      sp = round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
      data = VariableNodeData(zeros(dims, N),
                              zeros(dims,1), Int[], sp,
                              dims, false, 0, Int[], gt, softtype, false, false, dontmargin) #initialized
  end
  #
  setData!(v, data)
  # v.attributes["data"] = data

  # what about dedicated user variable metadata
  # p.factormetadata

  nothing
end

"""
    $(SIGNATURES)

Initialize a new Graphs.ExVertex which will be added to some factor graph.
"""
function addNewVarVertInGraph!(fgl::FactorGraph,
            vert::Graphs.ExVertex,
            id::Int,
            lbl::Symbol,
            ready::Int,
            smalldata  )
  #
  vert.attributes = Graphs.AttributeDict() #fg.v[fg.id]
  vert.attributes["label"] = string(lbl) #fg.v[fg.id]
  fgl.IDs[lbl] = id

  # used for cloudgraph solving
  vert.attributes["ready"] = ready
  vert.attributes["backendset"] = 0

  # store user metadata
  vert.attributes["smalldata"] = smalldata
  nothing
end



"""
$(SIGNATURES)

Add a variable node `lbl::Symbol` to `fg::FactorGraph`, as `softtype<:InferenceVariable`.

Example
-------

```julia
fg = initfg()
addVariable!(fg, :x0, Pose2)
```
"""
function addVariable!(fg::FactorGraph,
                      lbl::Symbol,
                      softtype::T;
                      N::Int=100,
                      autoinit::Bool=true,  # does init need to be separate from ready? TODO
                      ready::Int=1,
                      dontmargin::Bool=false,
                      labels::Vector{<:AbstractString}=String[],
                      api::DataLayerAPI=dlapi,
                      uid::Int=-1,
                      smalldata="",
                      checkduplicates::Bool=true  ) where {T <:InferenceVariable}
  #
  if checkduplicates
    if haskey(fg.IDs, lbl)
      @warn "Variable $lbl already exists in fg.sessionname=$(fg.sessionname).  Igoring this addVariable! call."
      return getVert(fg, lbl, api=api)
    end
  end

  currid = fg.id+1
  if uid==-1
    fg.id=currid
  else
    currid = uid
  end
  # dims = dims != -1 ? dims : st.dims
  lblstr = string(lbl)
  vert = ExVertex(currid,lblstr)
  addNewVarVertInGraph!(fg, vert, currid, lbl, ready, smalldata)
  # dlapi.setupvertgraph!(fg, vert, currid, lbl) #fg.v[currid]
  # dodims = fg.dimID+1
  setDefaultNodeData!(vert, 0, N, softtype.dims, initialized=!autoinit, softtype=softtype, dontmargin=dontmargin) # dodims

  vnlbls = union(string.(labels), softtype.labels, String["VARIABLE";])
  push!(vnlbls, fg.sessionname)

  api.addvertex!(fg, vert, labels=vnlbls)

  # fg.dimID+=dims # DONE -- drop this, rows indicate dimensions, move to last dimension
  push!(fg.nodeIDs, currid)

  # keep a fifo queue of incoming symbols
  push!(fg.fifo, lbl)

  vert
end

function addVariable!(fg::FactorGraph,
                      lbl::Symbol,
                      softtype::Type{<:InferenceVariable};
                      N::Int=100,
                      autoinit::Bool=true,
                      ready::Int=1,
                      dontmargin::Bool=false,
                      labels::Vector{<:AbstractString}=String[],
                      api::DataLayerAPI=dlapi,
                      uid::Int=-1,
                      # dims::Int=-1,
                      smalldata=""  )
  #
  sto = softtype()
  if :ut in fieldnames(typeof(sto))
    sto.ut != -9999999999 ? nothing : error("please define a microsecond time (;ut::Int64=___) for $(softtype)")
  end
  addVariable!(fg,
               lbl,
               sto,
               N=N,
               autoinit=autoinit,
               ready=ready,
               dontmargin=dontmargin,
               labels=labels,
               api=api,
               uid=uid,
               smalldata=smalldata  )
end



"""
    $(SIGNATURES)

Fetch the variable marginal sample points without the KDE bandwidth parameter.  Use getVertKDE to retrieve the full KDE object.
"""
function getVal(vA::Array{Graphs.ExVertex,1})
  @warn "getVal(::Vector{ExVertex}) is obsolete, use getVal.(ExVertex) instead."
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
  val = Array{Float64,2}(undef,rw, sc)
  for i in 1:(len-1)
      val[:,(cols[i]+1):cols[i+1]] = vals[i]
  end
  val[:,(cols[len]+1):cols[len+1]] = vals[len] # and the last one
  return val
end


"""
    $(SIGNATURES)

Prepare the particle arrays `ARR` to be used for approximate convolution.  This function ensures that ARR has te same dimensions among all the parameters.  Function returns with ARR[sfidx] pointing at newly allocated deepcopy of the existing values in getVal(Xi[.index==solvefor]).  Return values `sfidx` is the element in ARR where `Xi.index==solvefor` and `maxlen` is length of all (possibly resampled) `ARR` contained particles.  Note `Xi` is order sensitive.
"""
function prepareparamsarray!(ARR::Array{Array{Float64,2},1},
            Xi::Vector{Graphs.ExVertex},
            N::Int,
            solvefor::Int  )
  #
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

  # TODO --rather define reusable memory for the proposal
  # we are generating a proposal distribution, not direct replacement for existing memory and hence the deepcopy.
  if sfidx > 0 ARR[sfidx] = deepcopy(ARR[sfidx]) end
  return maxlen, sfidx
end

function parseusermultihypo(multihypo::Nothing)
  verts = Symbol[]
  mh = nothing
  return mh
end
function parseusermultihypo(multihypo::Union{Tuple,Vector{Float64}})
  mh = nothing
  if multihypo != nothing
    multihypo2 = Float64[multihypo...]
    # verts = Symbol.(multihypo[1,:])
    for i in 1:length(multihypo)
      if multihypo[i] > 0.999999
        multihypo2[i] = 0.0
      end
    end
    mh = Categorical(Float64[multihypo2...] )
  end
  return mh
end

# import IncrementalInference: prepgenericconvolution, convert

function prepgenericconvolution(
            Xi::Vector{Graphs.ExVertex},
            usrfnc::T;
            multihypo::Union{Nothing, Distributions.Categorical}=nothing,
            threadmodel=MultiThreaded  ) where {T <: FunctorInferenceType}
      # multiverts::Vector{Symbol}=Symbol[]
  #
  ARR = Array{Array{Float64,2},1}()
  maxlen, sfidx = prepareparamsarray!(ARR, Xi, 0, 0)
  fldnms = fieldnames(T) # typeof(usrfnc)
  zdim = T != GenericMarginal ? size(getSample(usrfnc, 2)[1],1) : 0
  certainhypo = multihypo != nothing ? collect(1:length(multihypo.p))[multihypo.p .== 0.0] : collect(1:length(Xi))
  ccw = CommonConvWrapper(
          usrfnc,
          zeros(1,0),
          zdim,
          ARR,
          specialzDim = sum(fldnms .== :zDim) >= 1,
          partial = sum(fldnms .== :partial) >= 1,
          hypotheses=multihypo,
          certainhypo=certainhypo,
          threadmodel=threadmodel
        )
  #
  for i in 1:Threads.nthreads()
    ccw.cpt[i].factormetadata.variableuserdata = []
    ccw.cpt[i].factormetadata.solvefor = :null
    for xi in Xi
      push!(ccw.cpt[i].factormetadata.variableuserdata, getData(xi).softtype)
    end
  end
  return ccw
end

function setDefaultFactorNode!(
      fgl::FactorGraph,
      vert::Graphs.ExVertex,
      Xi::Vector{Graphs.ExVertex},
      usrfnc::T;
      multihypo::Union{Nothing,Tuple,Vector{Float64}}=nothing,
      threadmodel=SingleThreaded  ) where {T <: Union{FunctorInferenceType, InferenceType}}
  #
  ftyp = typeof(usrfnc) # maybe this can be T
  # @show "setDefaultFactorNode!", usrfnc, ftyp, T
  mhcat = parseusermultihypo(multihypo)
  # gwpf = prepgenericwrapper(Xi, usrfnc, getSample, multihypo=mhcat)
  ccw = prepgenericconvolution(Xi, usrfnc, multihypo=mhcat, threadmodel=threadmodel)

  m = Symbol(ftyp.name.module)

  # experimental wip
  data_ccw = FunctionNodeData{CommonConvWrapper{T}}(Int[], false, false, Int[], m, ccw)
  vert.attributes["data"] = data_ccw

  # existing interface
  # data = FunctionNodeData{GenericWrapParam{T}}(Int[], false, false, Int[], m, gwpf)
  # vert.attributes["data"] = data

  nothing
end

function addNewFncVertInGraph!(fgl::FactorGraph, vert::Graphs.ExVertex, id::Int, lbl::Symbol, ready::Int)
  vert.attributes = Graphs.AttributeDict() #fg.v[fg.id]
  vert.attributes["label"] = lbl #fg.v[fg.id]
  # fgl.f[id] = vert #  -- not sure if this is required, using fg.g.vertices
  fgl.fIDs[lbl] = id # fg.id,

  # used for cloudgraph solving
  vert.attributes["ready"] = ready
  vert.attributes["backendset"] = 0

  # for graphviz drawing
  vert.attributes["shape"] = "point"
  vert.attributes["width"] = 0.2
  nothing
end
addNewFncVertInGraph!(fgl::FactorGraph, vert::Graphs.ExVertex, id::Int, lbl::T, ready::Int) where {T <: AbstractString} =
    addNewFncVertInGraph!(fgl,vert, id, Symbol(lbl), ready)

function isInitialized(vert::Graphs.ExVertex)::Bool
  return getData(vert).initialized
end
function isInitialized(fgl::FactorGraph, vsym::Symbol)::Bool
  # TODO, make cloudgraphs work and make faster by avoiding all the getVerts
  isInitialized(getVert(fgl, vsym))
end

"""
    $(SIGNATURES)

EXPERIMENTAL: initialize destination variable nodes based on this factor in factor graph, fg, generally called
during addFactor!. Destination factor is first (singletons) or second (dim 2 pairwise) variable vertex in Xi.
"""
function doautoinit!(fgl::FactorGraph,
                     Xi::Vector{Graphs.ExVertex};
                     api::DataLayerAPI=dlapi,
                     singles::Bool=true,
                     N::Int=100)
  # Mighty inefficient function, since we only need very select fields nearby from a few neighboring nodes
  # do double depth search for variable nodes
  # TODO this should maybe stay localapi only...
  for xi in Xi
    if !isInitialized(xi)
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
            if (isInitialized(vert2) && sum(useinitfct .== xifct) == 0 ) || length(xfneivarnodes) == 1
              # OR singleton  TODO get faster version of isInitialized for database version
              # println("adding $xifct to init factors list")
              push!(useinitfct, xifct)
            end
          end
        end
        # println("Consider all singleton (unary) factors to $vsym...")

        # calculate the predicted belief over $vsym
        pts = predictbelief(fgl, vsym, useinitfct, api=api)
        setValKDE!(xi, pts)
        getData(xi).initialized = true
        api.updatevertex!(fgl, xi, updateMAPest=false)
      end
    end
  end

  nothing
end

"""
    $(SIGNATURES)

Initialize destination variable nodes based on this factor in factor graph, fg, generally called
during addFactor!.  Destination factor is first (singletons) or second (dim 2 pairwise) variable vertex in Xi.
"""
function doautoinit!(fgl::FactorGraph,
                     xsyms::Vector{Symbol};
                     api::DataLayerAPI=dlapi,
                     singles::Bool=true,
                     N::Int=100)
  #
  verts = getVert.(fgl, xsyms, api=api)
  doautoinit!(fgl, verts, api=api, singles=singles, N=N)
end
function doautoinit!(fgl::FactorGraph,
                     xsym::Symbol;
                     api::DataLayerAPI=dlapi,
                     singles::Bool=true,
                     N::Int=100)
  #
  doautoinit!(fgl, [getVert(fgl, xsym, api=api);], api=api, singles=singles, N=N)
end

"""
    $(SIGNATURES)

Workaround function when first-version (factor graph based) auto initialization fails.  Usually occurs when using factors that have high connectivity to multiple variables.
"""
function manualinit!(fgl::FactorGraph, vert::Graphs.ExVertex, pX::BallTreeDensity)::Nothing
  setValKDE!(vert, pX)
  getData(vert).initialized = true
  # TODO must still update back to server is using api=dlapi
  nothing
end
function manualinit!(fgl::FactorGraph, sym::Symbol, pX::BallTreeDensity; api::DataLayerAPI=localapi)::Nothing
  vert = getVert(fgl, sym, api=api)
  manualinit!(fgl, vert, pX)
  nothing
end
function manualinit!(fgl::FactorGraph, sym::Symbol, usefcts::Vector{Symbol}; api::DataLayerAPI=localapi)::Nothing
  global localapi
  @warn "manual_init being used as a workaround for temporary autoinit issues."
  pts = predictbelief(fgl, sym, usefcts)
  vert = getVert(fgl, sym, api=api)
  Xpre = AMP.manikde!(pts, getSofttype(vert).manifolds )
  setValKDE!(vert, Xpre) # fgl, sym
  getData(fgl, sym).initialized = true
  nothing
end

function ensureAllInitialized!(fgl::FactorGraph; api::DataLayerAPI=dlapi)
  xx, xl = ls(fgl)
  allvarnodes = union(xx, xl)
  for vsym in allvarnodes
    if !isInitialized(fgl, vsym)
      @info "$vsym is not initialized, and will do so now..."
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

"""
    $(SIGNATURES)

Add factor with user defined type <: FunctorInferenceType to the factor graph object.  Define whether the automatic initialization of variables should be performed.  Use order sensitive `multihypo` keyword argument to define if any variables are related to data association uncertainty.
"""
function addFactor!(fgl::FactorGraph,
      Xi::Vector{Graphs.ExVertex},
      usrfnc::R;
      multihypo::Union{Nothing,Tuple,Vector{Float64}}=nothing,
      ready::Int=1,
      api::DataLayerAPI=dlapi,
      labels::Vector{T}=String[],
      uid::Int=-1,
      autoinit::Bool=true,
      threadmodel=SingleThreaded  ) where
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
  setDefaultFactorNode!(fgl, newvert, Xi, deepcopy(usrfnc), multihypo=multihypo, threadmodel=threadmodel)
  push!(fgl.factorIDs,currid)

  # TODO -- evaluate and streamline
  for vert in Xi
    push!(getData(newvert).fncargvID, vert.index)
    # push!(newvert.attributes["data"].fncargvID, vert.index)
  end

  fnlbls = deepcopy(labels)
  fnlbls = union(fnlbls, String["FACTOR";])
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
      multihypo::Union{Nothing,Tuple,Vector{Float64}}=nothing,
      ready::Int=1,
      api::DataLayerAPI=dlapi,
      labels::Vector{T}=String[],
      uid::Int=-1,
      autoinit::Bool=true,
      threadmodel=SingleThreaded  ) where
        {R <: Union{FunctorInferenceType, InferenceType},
         T <: AbstractString}
  #
  verts = Vector{Graphs.ExVertex}()
  for xi in xisyms
      push!( verts, api.getvertex(fgl,xi) )
  end
  addFactor!(fgl, verts, usrfnc, multihypo=multihypo, ready=ready, api=api, labels=labels, uid=uid, autoinit=autoinit, threadmodel=threadmodel )
end



"""
    $SIGNATURES

Delete factor and its edges.
"""
function deleteFactor!(fgl::FactorGraph, fsym::Symbol)
  fid = fgl.fIDs[fsym]
  eds = fgl.g.inclist[fid]
  alledsids = Int[]
  nedges = length(eds)
  for eds in fgl.g.inclist[fid]
    union!(alledsids, [eds.source.index; eds.target.index])
  end
  for edids in setdiff!(alledsids, fid)
    count = 0
    for eds in fgl.g.inclist[edids]
      count += 1
      if fid == eds.source.index || fid == eds.target.index
        deleteat!(fgl.g.inclist[edids], count)
        break
      end
    end
  end
  delete!(fgl.g.inclist, fid)
  fgl.g.nedges -= nedges
  delete!(fgl.g.vertices, fid)
  delete!(fgl.fIDs, fsym)
  deleteat!(fgl.factorIDs, findfirst(a -> a==fid, fgl.factorIDs))
  nothing
end

"""
    $SIGNATURES

Delete variables, and also the factors+edges if `andfactors=true` (default).
"""
function deleteVariable!(fgl::FactorGraph, vsym::Symbol; andfactors::Bool=true)
  vid = fgl.IDs[vsym]
  vert = fgl.g.vertices[vid]
  if andfactors
    for ne in Graphs.out_neighbors(vert, fgl.g)
      deleteFactor!(fgl, Symbol(ne.label))
    end
  end
  delete!(fgl.g.vertices, vid)
  delete!(fgl.IDs, vsym)
  deleteat!(fgl.nodeIDs, findfirst(a -> a==vid, fgl.nodeIDs))
  nothing
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
    permuteds = Vector{Int}(undef, lens)
    permutedsf = Vector{Int}(undef, lensf)
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
        p = cholfact(A'A,:U,Val(true))[:p] #,pivot=true
      elseif ordering==:qr
        q,r,p = qr(A, Val(true))
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
    # @show vert.label, getData(vert).BayesNetVertID
    if getData(vert).BayesNetVertID == 0
      fg.bnid+=1
      vert.attributes["data"].BayesNetVertID = p
      localapi.updatevertex!(fg, vert)
    else
      @warn "addBayesNetVerts -- something is very wrong, should not have a Bayes net vertex"
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
  # @info "adding marginal to"
  # for x in Xi
  #   @info "x.index=",x.index
  # end
  addFactor!(fg, Xi, genmarg, api=localapi, autoinit=false)
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
              @warn "removing edge $(edge.neo4jEdgeId), between $(m.index) and $(n.index)"
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
      @warn "removing vertex id=$(m.index)"
      localapi.deletevertex!(fgl,m)
    end
  end
  nothing
end

function buildBayesNet!(fg::FactorGraph, p::Array{Int,1})
    addBayesNetVerts!(fg, p)
    for v in p
      @info ""
      @info "Eliminating $(v)"
      @info "==============="
      @info ""
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
            # TODO -- validate !(sepNode.index in Si) vs. older !(sepNode in Si)
            if sepNode.index != v && !(sepNode.index in Si) # length(findin(sepNode.index, Si)) == 0
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

function getKDE(vnd::VariableNodeData)
  AMP.manikde!(getVal(vnd), getBW(vnd)[:,1], getSofttype(vnd).manifolds)
end


"""
    $(SIGNATURES)

Get KernelDensityEstimate kde estimate stored in variable node.
"""
function getKDE(v::Graphs.ExVertex)
  return getKDE(getData(v))
end

"""
    $(SIGNATURES)

Get KernelDensityEstimate kde estimate stored in variable node.
"""
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
function getKDE(fgl::FactorGraph, lbl::Symbol; api::DataLayerAPI=dlapi)
  return getVertKDE(fgl, lbl, api=api)
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

"""
    $(SIGNATURES)

Export a dot and pdf file drawn by Graphviz showing the factor graph.
"""
function writeGraphPdf(fgl::FactorGraph;
                       viewerapp::String="evince",
                       filepath::AS="/tmp/fg.pdf",
                       engine::AS="sfdp",
                       show::Bool=true ) where {AS <: AbstractString}
  #
  fgd = drawCopyFG(fgl)
  @info "Writing factor graph file"
  fext = split(filepath, '.')[end]
  fpwoext = split(filepath, '.')[end-1]
  dotfile = fpwoext*".dot"
  fid = open(dotfile,"w")
  write(fid,Graphs.to_dot(fgd.g))
  close(fid)
  show ? (@async run(`$(engine) $(dotfile) -T$(fext) -o $(filepath)`)) : nothing

  try
    viewerapp != nothing ? (@async run(`$(viewerapp) $(filepath)`)) : nothing
  catch e
    @warn "not able to show $(filepath) with viewerapp=$(viewerapp). Exception e=$(e)"
  end
  nothing
end


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


#
