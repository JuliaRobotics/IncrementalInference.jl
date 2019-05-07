# TODO: Remove this in time
import DistributedFactorGraphs.GraphsJl
const DFGGraphs = DistributedFactorGraphs.GraphsJl


reshapeVec2Mat(vec::Vector, rows::Int) = reshape(vec, rows, round(Int,length(vec)/rows))

# get vertex from factor graph according to label symbol "x1"
# getVert(dfg::FactorGraph, lbl::Symbol;, nt::Symbol=:var) = api.getvertex(fgl, lbl, nt=nt)
# Replacing with DFG.getVariable and DFG.getFactor

# """
#     $(TYPEDSIGNATURES)
#
# Get the `::Symbol` name for a node with `id::Int`.
# """
# getSym(fgl::FactorGraph, id::Int) = Symbol(fgl.g.vertices[id].label)

"""
    $SIGNATURES

Retrieve data structure stored in a node.
"""
getData(v::DFGFactor)::GenericFunctionNodeData = v.data
getData(v::DFGVariable; solveKey::Symbol=:default)::VariableNodeData = v.solverDataDict[solveKey]

"""
    $SIGNATURES

Retrieve data structure stored in a variable.
"""
function getVariableData(dfg::T, lbl::Symbol; solveKey::Symbol=:default)::VariableNodeData where {T <: AbstractDFG}
  return getData(getVariable(dfg, lbl, solveKey=solveKey))
end

"""
    $SIGNATURES

Retrieve data structure stored in a factor.
"""
function getFactorData(dfg::T, lbl::Symbol)::GenericFunctionNodeData where {T <: AbstractDFG}
  return getData(getFactor(dfg, lbl))
end
# TODO -- upgrade to dedicated memory location in Graphs.jl
# see JuliaArchive/Graphs.jl#233

# TODO: Intermediate for refactor. I'm sure we'll see this in 2024 though, it being 'temporary' and all :P
function setData!(v::DFGVariable, data::VariableNodeData; solveKey::Symbol=:default)::Nothing
  v.solverDataDict[solveKey] = data
  return nothing
end
function setData!(f::DFGFactor, data::GenericFunctionNodeData)::Nothing
  f.data = data
  return nothing
end

"""
   $(SIGNATURES)

Variable nodes softtype information holding a variety of meta data associated with the type of variable stored in that node of the factor graph.
"""
function getSofttype(vnd::VariableNodeData)
  return vnd.softtype
end
function getSofttype(v::DFGVariable; solveKey::Symbol=:default)
  return getSofttype(getData(v, solveKey=solveKey))
end


"""
    $SIGNATURES
Return the manifolds on which variable `sym::Symbol` is defined.
"""
getManifolds(v::DFGVariable; solveKey::Symbol=:default) = getSofttype(v, solveKey=solveKey).manifolds
function getManifolds(dfg::T, sym::Symbol; solveKey::Symbol=:default) where {T <: AbstractDFG}
  return getManifolds(getVariable(fgl, sym), solveKey=solveKey)
end

"""
    $(SIGNATURES)

Convenience function to get point values sampled i.i.d from marginal of `lbl` variable in the current factor graph.
"""
getVal(v::DFGVariable; solveKey::Symbol=:default) = v.solverDataDict[solveKey].val
getVal(v::DFGVariable, idx::Int; solveKey::Symbol=:default) = v.solverDataDict[solveKey].val[:,idx]
getVal(vnd::VariableNodeData) = vnd.val
getVal(vnd::VariableNodeData, idx::Int) = vnd.val[:, idx]
function getVal(dfg::T, lbl::Symbol; solveKey::Symbol=:default) where {T <: AbstractDFG}
  return getVariable(dfg, lbl).solverDataDict[solveKey]
end

"""
    $(SIGNATURES)

Get the number of points used for the current marginal belief estimate represtation for a particular variable in the factor graph.
"""
function getNumPts(v::DFGVariable; solveKey::Symbol=:default)::Tuple{Int64,Int64}
  return size(getData(v, solveKey=solveKey).val,2)
end

# TODO: Refactor - was is das?
function getfnctype(data::GenericFunctionNodeData)::Type
  if typeof(data).name.name == :VariableNodeData
    return VariableNodeData
  end
  return data.fnc.usrfnc!
end

function getfnctype(f::DFGFactor; solveKey::Symbol=:default)::Type
  data = getData(vertl, solveKey=solveKey)
  return getfnctype(data)
end

function getfnctype(dfg::T, lbl::Symbol; solveKey::Symbol=:default) where T <: AbstractDFG
  getfnctype(getFactor(dfg, exvertid, api=api))
end

function getBW(vnd::VariableNodeData)
  return vnd.bw
end

# setVal! assumes you will update values to database separate, this used for local graph mods only
function setVal!(v::DFGVariable, val::Array{Float64,2}; solveKey::Symbol=:default)::Nothing
  getData(v, solveKey=solveKey).val = val
  nothing
end
function getBWVal(v::DFGVariable; solveKey::Symbol=:default)
  return getData(v, solveKey=solveKey).bw
end
function setBW!(v::DFGVariable, bw::Array{Float64,2}; solveKey::Symbol=:default)::Nothing
  getData(v, solveKey=solveKey).bw = bw
  nothing
end
function setVal!(v::DFGVariable, val::Array{Float64,2}, bw::Array{Float64,2}; solveKey::Symbol=:default)::Nothing
  setVal!(v, val, solveKey=solveKey)
  setBW!(v, bw, solveKey=solveKey)
  nothing
end
function setVal!(v::DFGVariable, val::Array{Float64,2}, bw::Vector{Float64}; solveKey::Symbol=:default)
  setVal!(v, val, reshape(bw,length(bw),1), solveKey=solveKey)
  nothing
end
function setVal!(dfg::T, sym::Symbol, val::Array{Float64,2}; solveKey::Symbol=:default) where T <: AbstractDFG
  setVal!(getVariable(dfg, sym), val, solveKey=solveKey)
end

"""
    $SIGNATURES

Set the point centers and bandwidth parameters of a variable node, also set `isInitialized=true` if `setinit::Bool=true` (as per default).
"""
function setValKDE!(v::DFGVariable, val::Array{Float64,2}, setinit::Bool=true; solveKey::Symbol=:default)::Nothing
  # recover softtype information
  sty = getSofttype(v, solveKey=solveKey)
  # @show sty.manifolds
  #
  p = AMP.manikde!(val, sty.manifolds)
  setVal!(v,val,getBW(p)[:,1], solveKey=solveKey) # TODO -- this can be little faster
  setinit ? (getData(v, solveKey=solveKey).initialized = true) : nothing
  nothing
end
function setValKDE!(v::DFGVariable, em::EasyMessage, setinit::Bool=true; solveKey::Symbol=:default)::Nothing
  setVal!(v, em.pts, em.bws, solveKey=solveKey) # getBW(p)[:,1]
  setinit ? (getData(v, solveKey=solveKey).initialized = true) : nothing
  nothing
end
function setValKDE!(v::DFGVariable, p::BallTreeDensity, setinit::Bool=true; solveKey::Symbol=:default)
  pts = getPoints(p)
  setVal!(v, pts, getBW(p)[:,1], solveKey=solveKey) # BUG ...al!(., val, . ) ## TODO -- this can be little faster
  setinit ? (getData(v, solveKey==solveKey).initialized = true) : nothing
  nothing
end
function setValKDE!(dfg::T, sym::Symbol, p::BallTreeDensity; solveKey::Symbol=:default, setinit::Bool=true) where T <: AbstractDFG
  setValKDE!(getVert(dfg, sym), p, setinit, solveKey=solveKey)
  nothing
end
# TODO: Confirm this is supposed to be a variable?
setVal!(v::DFGVariable, em::EasyMessage; solveKey::Symbol=:default) = setValKDE!(v, em, solveKey=solveKey)
setVal!(v::DFGVariable, p::BallTreeDensity; solveKey::Symbol=:default) = setValKDE!(v, p, solveKey=solveKey)

"""
    $(SIGNATURES)

Construct a BallTreeDensity KDE object from an IIF.EasyMessage object.
"""
function kde!(em::EasyMessage)
  return AMP.manikde!(em.pts,em.bws, em.manifolds)
end


# TODO -- there should be a better way, without retrieving full vertex
function getOutNeighbors(dfg::T, v::V; needdata::Bool=false, ready::Union{Nothing, Int}=nothing, backendset::Union{Nothing, Int}=nothing)::Vector{Symbol} where {T <: AbstractDFG, V <: DFGNode}
  @warn "TODO: needdata is currently ignored. Symbols are returned."
  nodes = getNeighbors(dfg, v, ready=ready, backendset=backendset)
  return nodes
end
function getOutNeighbors(dfg::T, vertSym::Symbol; needdata::Bool=false, ready::Int=1, backendset::Int=1 )::Vector{Symbol} where {T <: AbstractDFG, V <: DFGNode}
  @warn "TODO: needdata is currently ignored. Symbols are returned."
  nodes = getNeighbors(dfg, vertSym, ready=ready, backendset=backendset)
  return nodes
end

# function updateFullVert!(fgl::FactorGraph, exvert::ExVertex;
#             api::DataLayerAPI=IncrementalInference.dlapi,
#             updateMAPest::Bool=false  )
#   #
#   @warn "use of updateFullVert! should be clarified for local or remote operations."
#   api.updatevertex!(fgl, exvert, updateMAPest=updateMAPest)
# end


function setDefaultNodeData!(v::DFGVariable,
                             dodims::Int,
                             N::Int,
                             dims::Int;
                             gt=nothing,
                             initialized::Bool=true,
                             dontmargin::Bool=false,
                             softtype=nothing)::Nothing
  # TODO review and refactor this function, exists as legacy from pre-v0.3.0
  # this should be the only function allocating memory for the node points (unless number of points are changed)
  data = nothing
  if initialized
    # if size(initval,2) < N && size(initval, 1) == dims
    #   @warn "setDefaultNodeData! -- deprecated use of stdev."
    #   p = AMP.manikde!(initval,diag(stdev), softtype.manifolds);
    #   pN = resample(p,N)
    # if size(initval,2) < N && size(initval, 1) != dims
      # @info "Node value memory allocated but not initialized"
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
    setSolverData(v, VariableNodeData(pNpts,
                            gbw2, Int[], sp,
                            dims, false, 0, Int[], gt, softtype, true, false, dontmargin))
  else
    sp = round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
    setSolverData(v, VariableNodeData(zeros(dims, N),
                            zeros(dims,1), Int[], sp,
                            dims, false, 0, Int[], gt, softtype, false, false, dontmargin))
  end
  return nothing
end

# """
#     $(SIGNATURES)
#
# Initialize a new Graphs.ExVertex which will be added to some factor graph.
# """
# function addNewVarVertInGraph!(fgl::FactorGraph,
#             vert::Graphs.ExVertex,
#             id::Int,
#             lbl::Symbol,
#             ready::Int,
#             smalldata  )
#   @error("NEIN NEIN ES IS VERBOTEN! DEPRECATED JA")
#   #
#   # vert.attributes = Graphs.AttributeDict() #fg.v[fg.id]
#   # vert.attributes["label"] = string(lbl) #fg.v[fg.id]
#   # fgl.IDs[lbl] = id
#   #
#   # # used for cloudgraph solving
#   # vert.attributes["ready"] = ready
#   # vert.attributes["backendset"] = 0
#   #
#   # # store user metadata
#   # vert.attributes["smalldata"] = smalldata
#   # nothing
# end



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
function addVariable!(dfg::G,
                      lbl::Symbol,
                      softtype::T;
                      N::Int=100,
                      autoinit::Bool=true,  # does init need to be separate from ready? TODO
                      ready::Int=1,
                      dontmargin::Bool=false,
                      labels::Vector{Symbol}=Symbol[],
                      smalldata="",
                      checkduplicates::Bool=true  )::DFGVariable where
                      {G <: AbstractDFG,
                       T <: InferenceVariable}
  #
  # if checkduplicates
  #   if haskey(fg.IDs, lbl)
  #     @warn "Variable $lbl already exists in fg.sessionname=$(fg.sessionname).  Igoring this addVariable! call."
  #     return getVert(fg, lbl, api=api)
  #   end
  # end

  # currid = fg.id+1
  # if uid==-1
  #   fg.id=currid
  # else
  #   currid = uid
  # end
  v = DFGVariable(lbl)
  v.ready = ready
  # v.backendset = backendset
  v.tags = union(labels, Symbol.(softtype.labels), [:VARIABLE])
  v.smallData = smalldata
  # dims = dims != -1 ? dims : st.dims
  # lblstr = string(lbl)
  # vert = ExVertex(currid,lblstr)
  # addNewVarVertInGraph!(fg, vert, currid, lbl, ready, smalldata)
  # dlapi.setupvertgraph!(fg, vert, currid, lbl) #fg.v[currid]
  # dodims = fg.dimID+1
  setDefaultNodeData!(v, 0, N, softtype.dims, initialized=!autoinit, softtype=softtype, dontmargin=dontmargin) # dodims
  # TODO: Make generic
  DFGGraphs.addVariable!(dfg, v)

  # vnlbls = union(string.(labels), softtype.labels, String["VARIABLE";])
  # push!(vnlbls, fg.sessionname)

  # api.addvertex!(fg, vert, labels=vnlbls)
  #
  # # fg.dimID+=dims # DONE -- drop this, rows indicate dimensions, move to last dimension
  # push!(fg.nodeIDs, currid)

  # keep a fifo queue of incoming symbols
  # NOTE: Now using getAddHistory(dfg) for this - it's done internally
  # push!(fg.fifo, lbl)

  return v
end


function addVariable!(dfg::G,
                      lbl::Symbol,
                      softtype::Type{<:InferenceVariable};
                      N::Int=100,
                      autoinit::Bool=true,
                      ready::Int=1,
                      dontmargin::Bool=false,
                      labels::Vector{Symbol}=Symbol[],
                      smalldata="")::DFGVariable where
                      {G <: AbstractDFG} #
  sto = softtype()
  #TODO: Refactor
  if :ut in fieldnames(typeof(sto))
    sto.ut != -9999999999 ? nothing : error("please define a microsecond time (;ut::Int64=___) for $(softtype)")
  end
  return addVariable!(dfg,
               lbl,
               sto,
               N=N,
               autoinit=autoinit,
               ready=ready,
               dontmargin=dontmargin,
               labels=labels,
               smalldata=smalldata  )
end



"""
    $(SIGNATURES)

Fetch the variable marginal sample points without the KDE bandwidth parameter.  Use getVertKDE to retrieve the full KDE object.
"""
function getVal(vA::Vector{DFGVariable}, solveKey::Symbol=:default)::Array{Float64, 2}
  @warn "getVal(::Vector{ExVertex}) is obsolete, use getVal.(ExVertex) instead."
  len = length(vA)
  vals = Array{Array{Float64,2},1}()
  cols = Array{Int,1}()
  push!(cols,0)
  rows = Array{Int,1}()
  for v in vA
      push!(vals, getVal(v, solveKey=solveKey))
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

Prepare the particle arrays `ARR` to be used for approximate convolution.
This function ensures that ARR has te same dimensions among all the parameters.
Function returns with ARR[sfidx] pointing at newly allocated deepcopy of the
existing values in getVal(Xi[.label==solvefor]).
Return values `sfidx` is the element in ARR where `Xi.label==solvefor` and
`maxlen` is length of all (possibly resampled) `ARR` contained particles.
Note `Xi` is order sensitive.
Note for initialization, solveFor = Nothing.
"""
function prepareparamsarray!(ARR::Array{Array{Float64,2},1},
            Xi::Vector{DFGVariable},
            N::Int,
            solvefor::Union{Nothing, Symbol}  ) # TODO: Confirm we can use symbols here
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
    if xi.label == solvefor
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
            Xi::Vector{DFGVariable},
            usrfnc::T;
            multihypo::Union{Nothing, Distributions.Categorical}=nothing,
            threadmodel=MultiThreaded  ) where {T <: FunctorInferenceType}
      # multiverts::Vector{Symbol}=Symbol[]
  #
  ARR = Array{Array{Float64,2},1}()
  maxlen, sfidx = prepareparamsarray!(ARR, Xi, 0, nothing) # Nothing for init.
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
      dfg::G,
      factor::DFGFactor,
      Xi::Vector{DFGVariable},
      usrfnc::T;
      multihypo::Union{Nothing,Tuple,Vector{Float64}}=nothing,
      threadmodel=SingleThreaded  )::GenericFunctionNodeData where
      {G <: AbstractDFG, T <: Union{FunctorInferenceType, InferenceType}}
  #
  ftyp = T # typeof(usrfnc) # maybe this can be T
  # @show "setDefaultFactorNode!", usrfnc, ftyp, T
  mhcat = parseusermultihypo(multihypo)
  # gwpf = prepgenericwrapper(Xi, usrfnc, getSample, multihypo=mhcat)
  ccw = prepgenericconvolution(Xi, usrfnc, multihypo=mhcat, threadmodel=threadmodel)

  # experimental wip
  data_ccw = FunctionNodeData{CommonConvWrapper{T}}(Int[], false, false, Int[], Symbol(ftyp.name.module), ccw)
  factor.data = data_ccw

  # existing interface
  # data = FunctionNodeData{GenericWrapParam{T}}(Int[], false, false, Int[], m, gwpf)
  # vert.attributes["data"] = data

  return factor.data
end

# """
#     $SIGNATURES
#
# Return the current maximum vertex ID number in the factor graph fragment `fgl`.
# """
# function getMaxVertId(fgl::FactorGraph)
#   ids = collect(keys(fgl.g.vertices))
#   return maximum(ids)
# end

# function addNewFncVertInGraph!(fgl::FactorGraph, vert::Graphs.ExVertex, id::Int, lbl::Symbol, ready::Int)
#   vert.attributes = Graphs.AttributeDict() #fg.v[fg.id]
#   vert.attributes["label"] = lbl #fg.v[fg.id]
#   # fgl.f[id] = vert #  -- not sure if this is required, using fg.g.vertices
#   fgl.fIDs[lbl] = id # fg.id,
#
#   # used for cloudgraph solving
#   vert.attributes["ready"] = ready
#   vert.attributes["backendset"] = 0
#
#   # for graphviz drawing
#   vert.attributes["shape"] = "point"
#   vert.attributes["width"] = 0.2
#   nothing
# end
# addNewFncVertInGraph!(fgl::FactorGraph, vert::Graphs.ExVertex, id::Int, lbl::T, ready::Int) where {T <: AbstractString} =
#     addNewFncVertInGraph!(fgl,vert, id, Symbol(lbl), ready)

"""
    $SIGNATURES

Returns state of vertex data `.initialized` flag.

Notes:
- used by both factor graph variable and Bayes tree clique logic.
TODO: Refactor
"""
function isInitialized(vert::Graphs.ExVertex)::Bool
  return getData(vert).initialized
end
function isInitialized(vert::DFGVariable)::Bool
  return getData(vert).initialized
end
function isInitialized(dfg::T, vsym::Symbol)::Bool where T <: AbstractDFG
  return isInitialized(getVariable(dfg, vsym))
end

"""
    $SIGNATURES

Return `(::Bool, ::OKVarlist, ::NotOkayVarList)` on whether all other variables (besides `loovar::Symbol`)
attached to factor `fct::Symbol` are all initialized -- i.e. `fct` is usable.

Development Notes
* TODO get faster version of isInitialized for database version
"""
function factorCanInitFromOtherVars(dfg::T,
                                    fct::Symbol,
                                    loovar::Symbol)::Tuple{Bool, Vector{Symbol}, Vector{Symbol}} where T <: AbstractDFG
  #
  # all variables attached to this factor
  varsyms = getNeighbors(dfg, fct)

  # list of factors to use in init operation
  useinitfct = Symbol[]
  faillist = Symbol[]
  for vsym in varsyms
    xi = DFGGraphs.getVariable(dfg, vsym)
    if (isInitialized(xi) && sum(useinitfct .== fct) == 0 ) || length(varsyms) == 1
      push!(useinitfct, fct)
    end
  end

  return (length(useinitfct)==length(varsyms)&&length(faillist)==0,
          useinitfct,
          faillist   )
end

"""
    $(SIGNATURES)

EXPERIMENTAL: initialize target variable `xi` based on connected factors in the
factor graph `fgl`.  Possibly called from `addFactor!`, or `doCliqAutoInitUp!`.

Development Notes:
> Target factor is first (singletons) or second (dim 2 pairwise) variable vertex in `xi`.
* TODO use DFG properly with local operations and DB update at end.
* TODO get faster version of `isInitialized` for database version.
* TODO: Persist this back if we want to here.
"""
function doautoinit!(dfg::T,
                     xi::DFGVariable;
                     singles::Bool=true,
                     N::Int=100)::Bool where T <: AbstractDFG
  #
  didinit = false
  # don't initialize a variable more than once
  if !isInitialized(xi)
    # get factors attached to this variable xi
    vsym = Symbol(xi.label)
    neinodes = DFGGraphs.getNeighbors(dfg, vsym)
    # proceed if has more than one neighbor OR even if single factor
    if (length(neinodes) > 1 || singles) # && !isInitialized(xi)
      # Which of the factors can be used for initialization
      useinitfct = Symbol[]
      # Consider factors connected to $vsym...
      for xifct in neinodes
        canuse, usefct, notusevars = factorCanInitFromOtherVars(dfg, xifct, vsym)
        useinitfct = union(useinitfct, usefct)
      end
      # println("Consider all singleton (unary) factors to $vsym...")
      # calculate the predicted belief over $vsym
      if length(useinitfct) > 0
        pts = predictbelief(dfg, vsym, useinitfct)
        setValKDE!(xi, pts)
        getData(xi).initialized = true
        # TODO: Persist this back if we want to here.
        # api.updatevertex!(dfg, xi, updateMAPest=false)
        didinit = true
      end
    end
  end
  return didinit
end

function doautoinit!(dfg::T,
                     Xi::Vector{DFGVariable};
                     singles::Bool=true,
                     N::Int=100)::Bool where T <: AbstractDFG
  #
  #
  # Mighty inefficient function, since we only need very select fields nearby from a few neighboring nodes
  # do double depth search for variable nodes

  didinit = true

  # loop over all requested variables that must be initialized
  for xi in Xi
    didinit &= doautoinit!(dfg, xi, singles=singles, N=N)
  end
  return didinit
end

function doautoinit!(dfg::T,
                     xsyms::Vector{Symbol};
                     singles::Bool=true,
                     N::Int=100)::Bool where T <: AbstractDFG
  #
  verts = getVariable.(dfg, xsyms)
  return doautoinit!(dfg, verts, singles=singles, N=N)
end
function doautoinit!(dfg::T,
                     xsym::Symbol;
                     singles::Bool=true,
                     N::Int=100)::Bool where T <: AbstractDFG
  #
  return doautoinit!(dfg, [getVariable(dfg, xsym);], singles=singles, N=N)
end

"""
    $(SIGNATURES)

Workaround function when first-version (factor graph based) auto initialization fails.  Usually occurs when using factors that have high connectivity to multiple variables.
"""
function manualinit!(dfg::T, vert::DFGVariable, pX::BallTreeDensity)::Nothing where T <: AbstractDFG
  setValKDE!(vert, pX)
  getData(vert).initialized = true
  return nothing
end
function manualinit!(dfg::T, sym::Symbol, pX::BallTreeDensity)::Nothing where T <: AbstractDFG
  vert = getVariable(dfg, sym)
  manualinit!(dfg, vert, pX)
  return nothing
end
function manualinit!(dfg::T, sym::Symbol, usefcts::Vector{Symbol})::Nothing where T <: AbstractDFG
  @warn "manual_init being used as a workaround for temporary autoinit issues."
  pts = predictbelief(dfg, sym, usefcts)
  vert = getVert(dfg, sym, api=api)
  Xpre = AMP.manikde!(pts, getSofttype(vert).manifolds )
  setValKDE!(vert, Xpre) # dfg, sym
  getData(dfg, sym).initialized = true
  return nothing
end

function ensureAllInitialized!(dfg::T) where T <: AbstractDFG
  allvarnodes = getVariables(dfg)
  for var in allvarnodes
    if !isInitialized(var)
      @info "$(var.label) is not initialized, and will do so now..."
      doautoinit!(dfg, [var;], singles=true)
    end
  end
  nothing
end

function assembleFactorName(dfg::T, Xi::Vector{DFGVariable}; maxparallel::Int=50)::String where T <: AbstractDFG
  existingFactorLabels = getFactorIds(dfg)
  existingFactorLabelDict = Dict(existingFactorLabels .=> existingFactorLabels)
  namestring = ""
  for vert in Xi #f.Xi
    namestring = string(namestring,vert.label)
  end
  for i in 1:maxparallel
    tempnm = string(namestring, "f$i")
    if !haskey(existingFactorLabelDict, Symbol(tempnm))
      namestring = tempnm
      break
    end
    i != maxparallel ? nothing : error("Cannot currently add more than $(maxparallel) factors in parallel, please open an issue if this is too restrictive.")
  end
  return namestring
end

"""
    $(SIGNATURES)

Add factor with user defined type <: FunctorInferenceType to the factor graph
object. Define whether the automatic initialization of variables should be
performed.  Use order sensitive `multihypo` keyword argument to define if any
variables are related to data association uncertainty.
"""
function addFactor!(
      dfg::G,
      Xi::Vector{DFGVariable},
      usrfnc::R;
      multihypo::Union{Nothing,Tuple,Vector{Float64}}=nothing,
      ready::Int=1,
      labels::Vector{Symbol}=Symbol[],
      autoinit::Bool=true,
      threadmodel=SingleThreaded  ) where
        {G <: AbstractDFG,
         R <: Union{FunctorInferenceType, InferenceType}}
  #
  namestring = assembleFactorName(dfg, Xi)
  newFactor = DFGFactor{CommonConvWrapper{R}, Symbol}(Symbol(namestring))
  newFactor.tags = union(labels, [:FACTOR]) # TODO: And session info
  # addNewFncVertInGraph!(fgl, newvert, currid, namestring, ready)
  newData = setDefaultFactorNode!(dfg, newFactor, Xi, deepcopy(usrfnc), multihypo=multihypo, threadmodel=threadmodel)

  # TODO: Need to remove this...
  for vert in Xi
    push!(newData.fncargvID, vert._internalId) # YUCK :/
  end

  success = DFGGraphs.addFactor!(dfg, Xi, newFactor)

  # TODO: change this operation to update a conditioning variable
  autoinit && doautoinit!(dfg, Xi, singles=false)

  return newFactor
end
function addFactor!(
      dfg::G,
      xisyms::Vector{Symbol},
      usrfnc::R;
      multihypo::Union{Nothing,Tuple,Vector{Float64}}=nothing,
      ready::Int=1,
      labels::Vector{Symbol}=Symbol[],
      autoinit::Bool=true,
      threadmodel=SingleThreaded  ) where
        {G <: AbstractDFG,
         R <: Union{FunctorInferenceType, InferenceType}}
  #
  verts = map(vid -> DFGGraphs.getVariable(dfg, vid), xisyms)
  addFactor!(dfg, verts, usrfnc, multihypo=multihypo, ready=ready, labels=labels, autoinit=autoinit, threadmodel=threadmodel )
end



# """
#     $SIGNATURES
#
# Delete factor and its edges.
# """
# function deleteFactor!(fgl::FactorGraph, fsym::Symbol)
#   fid = fgl.fIDs[fsym]
#   eds = fgl.g.inclist[fid]
#   alledsids = Int[]
#   nedges = length(eds)
#   for eds in fgl.g.inclist[fid]
#     union!(alledsids, [eds.source.index; eds.target.index])
#   end
#   for edids in setdiff!(alledsids, fid)
#     count = 0
#     for eds in fgl.g.inclist[edids]
#       count += 1
#       if fid == eds.source.index || fid == eds.target.index
#         deleteat!(fgl.g.inclist[edids], count)
#         break
#       end
#     end
#   end
#   delete!(fgl.g.inclist, fid)
#   fgl.g.nedges -= nedges
#   delete!(fgl.g.vertices, fid)
#   delete!(fgl.fIDs, fsym)
#   deleteat!(fgl.factorIDs, findfirst(a -> a==fid, fgl.factorIDs))
#   nothing
# end

# """
#     $SIGNATURES
#
# Delete variables, and also the factors+edges if `andfactors=true` (default).
# """
# function deleteVariable!(fgl::FactorGraph, vsym::Symbol; andfactors::Bool=true)
#   vid = fgl.IDs[vsym]
#   vert = fgl.g.vertices[vid]
#   if andfactors
#     for ne in Graphs.out_neighbors(vert, fgl.g)
#       deleteFactor!(fgl, Symbol(ne.label))
#     end
#   end
#   delete!(fgl.g.vertices, vid)
#   delete!(fgl.IDs, vsym)
#   deleteat!(fgl.nodeIDs, findfirst(a -> a==vid, fgl.nodeIDs))
#   nothing
# end


function prtslperr(s)
  println(s)
  sleep(0.1)
  error(s)
end

### TODO: TO BE REFACTORED FOR DFG

# for computing the Bayes Net-----------------------------------------------------
function getEliminationOrder(dfg::G; ordering::Symbol=:qr) where G <: AbstractDFG
  s = getVariableIds(dfg)
  lens = length(s)
  sf = getFactorIds(dfg)
  lensf = length(sf)
  # adjm, dictpermu = adjacency_matrix(fg.g,returnpermutation=true)
  @show adjm = getAdjacencyMatrix(dfg)
  @error "Getting adjacency matrix here and it's new behavior, this is probably going to break!"
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
  return dictpermu[permuteds[p]] # fg.nodeIDs[p]
end


# lets create all the vertices first and then deal with the elimination variables thereafter
function addBayesNetVerts!(dfg::G, elimOrder::Array{Int,1}) where G <: AbstractDFG
  @error "This method works on ID's, we should use labels, needs a refactor."
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
      getData(tuv).eliminated = true # fg.v[v].
      localapi.updatevertex!(fg, tuv)

      # TODO -- remove links from current vertex to any marginals
      rmVarFromMarg(fg, vert, gm)

      #add marginal on remaining variables... ? f(xyz) = f(x | yz) f(yz)
      # new function between all Si
      addChainRuleMarginal!(fg, Si)

    end
    nothing
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
function writeGraphPdf(dfg::G;
                       viewerapp::String="evince",
                       filepath::AS="/tmp/fg.pdf",
                       engine::AS="sfdp",
                       show::Bool=true ) where {G <: AbstractDFG, AS <: AbstractString}
  #
  # fgd = drawCopyFG(fgl)
  @info "Writing factor graph file"
  fext = split(filepath, '.')[end]
  fpwoext = split(filepath, '.')[end-1]
  dotfile = fpwoext*".dot"
  # fid = open(dotfile,"w")
  # write(fid,Graphs.to_dot(fgd.g))
  # close(fid)
  toDotFile(dfg, dotfile)
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
