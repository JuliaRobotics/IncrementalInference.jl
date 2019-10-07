reshapeVec2Mat(vec::Vector, rows::Int) = reshape(vec, rows, round(Int,length(vec)/rows))


import DistributedFactorGraphs: getData
"""
    $SIGNATURES

Retrieve data structure stored in a node.
"""
getData(v::DFGFactor)::GenericFunctionNodeData = v.data
getData(v::DFGVariable; solveKey::Symbol=:default)::VariableNodeData = v.solverDataDict[solveKey]
# For Bayes tree
getData(v::Graphs.ExVertex) = v.attributes["data"]
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
# For Bayes tree
function setData!(v::Graphs.ExVertex, data)
  # this is a memory gulp without replacement, old attr["data"] object is left to gc
  v.attributes["data"] = data
  nothing
end

## has been moved to DFG
import DistributedFactorGraphs: getSofttype
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
getManifolds(vd::VariableNodeData) = getSofttype(vd).manifolds
getManifolds(v::DFGVariable; solveKey::Symbol=:default) = getManifolds(getData(v, solveKey=solveKey))
function getManifolds(dfg::G, sym::Symbol; solveKey::Symbol=:default) where {G <: AbstractDFG}
  return getManifolds(getVariable(dfg, sym), solveKey=solveKey)
end
getManifolds(vartype::InferenceVariable) = vartype.manifolds
getManifolds(vartype::Type{<: InferenceVariable}) = getManifolds(vartype())

"""
    $(SIGNATURES)

Convenience function to get point values sampled i.i.d from marginal of `lbl` variable in the current factor graph.
"""
getVal(v::DFGVariable; solveKey::Symbol=:default) = v.solverDataDict[solveKey].val
getVal(v::DFGVariable, idx::Int; solveKey::Symbol=:default) = v.solverDataDict[solveKey].val[:,idx]
getVal(vnd::VariableNodeData) = vnd.val
getVal(vnd::VariableNodeData, idx::Int) = vnd.val[:, idx]
function getVal(dfg::T, lbl::Symbol; solveKey::Symbol=:default) where {T <: AbstractDFG}
  return getVariable(dfg, lbl).solverDataDict[solveKey].val
end

"""
    $(SIGNATURES)

Get the number of points used for the current marginal belief estimate represtation for a particular variable in the factor graph.
"""
function getNumPts(v::DFGVariable; solveKey::Symbol=:default)::Int
  return size(getData(v, solveKey=solveKey).val,2)
end

import DistributedFactorGraphs: getfnctype
# TODO: Refactor - was is das?
function getfnctype(data::GenericFunctionNodeData)
  if typeof(data).name.name == :VariableNodeData
    return VariableNodeData
  end
  return data.fnc.usrfnc!
end

function getfnctype(fact::DFGFactor; solveKey::Symbol=:default)
  data = getData(fact) # TODO , solveKey=solveKey)
  return getfnctype(data)
end

function getfnctype(dfg::T, lbl::Symbol; solveKey::Symbol=:default) where T <: AbstractDFG
  getfnctype(getFactor(dfg, exvertid))
end

function getBW(vnd::VariableNodeData)
  return vnd.bw
end

# setVal! assumes you will update values to database separate, this used for local graph mods only
function getBWVal(v::DFGVariable; solveKey::Symbol=:default)
  return getData(v, solveKey=solveKey).bw
end
function setBW!(vd::VariableNodeData, bw::Array{Float64,2}; solveKey::Symbol=:default)::Nothing
  vd.bw = bw
  nothing
end
function setBW!(v::DFGVariable, bw::Array{Float64,2}; solveKey::Symbol=:default)::Nothing
  setBW!(getData(v, solveKey=solveKey), bw)
  nothing
end

function setVal!(vd::VariableNodeData, val::Array{Float64,2})::Nothing
    vd.val = val
    nothing
end
function setVal!(v::DFGVariable, val::Array{Float64,2}; solveKey::Symbol=:default)::Nothing
    setVal!(getData(v, solveKey=solveKey), val)
    nothing
end
function setVal!(vd::VariableNodeData, val::Array{Float64,2}, bw::Array{Float64,2})::Nothing
    setVal!(vd, val)
    setBW!(vd, bw)
    nothing
end
function setVal!(v::DFGVariable, val::Array{Float64,2}, bw::Array{Float64,2}; solveKey::Symbol=:default)::Nothing
  setVal!(v, val, solveKey=solveKey)
  setBW!(v, bw, solveKey=solveKey)
  nothing
end
function setVal!(vd::VariableNodeData, val::Array{Float64,2}, bw::Vector{Float64}; solveKey::Symbol=:default)
  setVal!(vd, val, reshape(bw,length(bw),1))
  nothing
end
function setVal!(v::DFGVariable, val::Array{Float64,2}, bw::Vector{Float64}; solveKey::Symbol=:default)
  setVal!(getData(v, solveKey=solveKey), val, bw)
  nothing
end
function setVal!(dfg::T, sym::Symbol, val::Array{Float64,2}; solveKey::Symbol=:default) where T <: AbstractDFG
  setVal!(getVariable(dfg, sym), val, solveKey=solveKey)
end

"""
    $SIGNATURES

Set the point centers and bandwidth parameters of a variable node, also set `isInitialized=true` if `setinit::Bool=true` (as per default).

Notes
- `initialized` is used for initial solve of factor graph where variables are not yet initialized.
- `inferdim` is used to identify if the initialized was only partial.
"""
function setValKDE!(vd::VariableNodeData,
                    pts::Array{Float64,2},
                    bws::Vector{Float64},
                    setinit::Bool=true,
                    inferdim::Float64=0 )::Nothing
  #
  setVal!(vd, pts, bws) # BUG ...al!(., val, . ) ## TODO -- this can be a little faster
  setinit ? (vd.initialized = true) : nothing
  vd.inferdim = inferdim
  nothing
end

function setValKDE!(vd::VariableNodeData,
                    p::BallTreeDensity,
                    setinit::Bool=true,
                    inferdim::Union{Float32, Float64, Int32, Int64}=0 )
  #
  pts = getPoints(p)
  bws = getBW(p)[:,1]
  setValKDE!(vd,pts,bws,setinit,inferdim )
  nothing
end

function setValKDE!(vd::VariableNodeData,
                    val::Array{Float64,2},
                    setinit::Bool=true,
                    inferdim::Float64=0 )::Nothing
  # recover softtype information
  sty = getSofttype(vd)
  p = AMP.manikde!(val, sty.manifolds)
  setValKDE!(vd, p, setinit, inferdim)
  # setVal!(vd,val,getBW(p)[:,1]) # TODO -- this can be little faster
  # setinit ? (vd.initialized = true) : nothing
  # vd.inferdim = inferdim
  nothing
end

function setValKDE!(v::DFGVariable,
                    val::Array{Float64,2},
                    bws::Vector{Float64},
                    setinit::Bool=true,
                    inferdim::Float64=0;
                    solveKey::Symbol=:default)::Nothing
  # recover softtype information
  setValKDE!(getData(v, solveKey=solveKey),val, bws, setinit, inferdim )

  nothing
end

function setValKDE!(v::DFGVariable,
                    val::Array{Float64,2},
                    setinit::Bool=true,
                    inferdim::Float64=0;
                    solveKey::Symbol=:default)::Nothing
  # recover softtype information
  setValKDE!(getData(v, solveKey=solveKey),val, setinit, inferdim )
  # sty = getSofttype(v, solveKey=solveKey)
  # p = AMP.manikde!(val, sty.manifolds)
  # setVal!(v,val,getBW(p)[:,1], solveKey=solveKey) # TODO -- this can be little faster
  # setinit ? (getData(v, solveKey=solveKey).initialized = true) : nothing
  # getData(v).inferdim = inferdim
  nothing
end
function setValKDE!(v::DFGVariable,
                    em::EasyMessage,
                    setinit::Bool=true,
                    # inferdim::Union{Float32, Float64, Int32, Int64}=0;
                    solveKey::Symbol=:default  )::Nothing
  #
  setValKDE!(v, em.pts, em.bws, setinit, em.inferdim, solveKey=solveKey) # getBW(p)[:,1]
  # setinit ? (getData(v, solveKey=solveKey).initialized = true) : nothing
  # getData(v).inferdim = inferdim
  nothing
end
function setValKDE!(v::DFGVariable,
                    p::BallTreeDensity,
                    setinit::Bool=true,
                    inferdim::Union{Float32, Float64, Int32, Int64}=0;
                    solveKey::Symbol=:default  )
  #
  setValKDE!(getData(v,solveKey=solveKey),p,setinit,Float64(inferdim))
  # pts = getPoints(p)
  # setVal!(v, pts, getBW(p)[:,1], solveKey=solveKey) # BUG ...al!(., val, . ) ## TODO -- this can be little faster
  # setinit ? (getData(v, solveKey=solveKey).initialized = true) : nothing
  # getData(v).inferdim = inferdim
  nothing
end
function setValKDE!(dfg::G,
                    sym::Symbol,
                    p::BallTreeDensity,
                    setinit::Bool=true,
                    inferdim::Union{Float32, Float64, Int32, Int64}=0;
                    solveKey::Symbol=:default  ) where G <: AbstractDFG
    #
    setValKDE!(getVariable(dfg, sym), p, setinit, inferdim, solveKey=solveKey)
    nothing
end

# TODO: Confirm this is supposed to be a variable?
function setVal!(v::DFGVariable, em::EasyMessage; solveKey::Symbol=:default)
    @warn "setVal! deprecated, use setValKDE! instead"
    setValKDE!(v, em, solveKey=solveKey)
end
function setVal!(v::DFGVariable, p::BallTreeDensity; solveKey::Symbol=:default)
    @warn "setVal! deprecated, use setValKDE! instead"
    setValKDE!(v, p, solveKey=solveKey)
end

"""
    $(SIGNATURES)

Construct a BallTreeDensity KDE object from an IIF.EasyMessage object.

Related

manikde!, getKDE, getKDEMax, getKDEMean, EasyMessage
"""
function kde!(em::EasyMessage)
  return AMP.manikde!(em.pts, em.bws, em.manifolds)
end



"""
    $SIGNATURES

Set variable initialized status.
"""
function setVariableInitialized!(varid::VariableNodeData,
                                 status::Bool)
  #
  varid.initialized = status
end

setVariableInitialized!(vari::DFGVariable, status::Bool) = setVariableInitialized!(getData(vari), status)


"""
    $SIGNATURES

Set method for the inferred dimension value in a variable.
"""
setVariableInferDim!(varid::VariableNodeData, val::Real) = varid.inferdim = convert(Float64,val)
setVariableInferDim!(vari::DFGVariable, val::Real) = setVariableInferDim!(getData(vari), val)


"""
    $SIGNATURES

Reset the solve state of a variable to uninitialized/unsolved state.
"""
function resetVariable!(varid::VariableNodeData;
                        solveKey::Symbol=:default  )::Nothing
  #
  val = getKDE(varid)
  pts = getPoints(val)
  # TODO not all manifolds will initialize to zero
  fill!(pts, 0.0)
  pn = manikde!(pts, zeros(KDE.Ndim(val)), getManifolds(varid))
  setValKDE!(varid, pn, false, 0.0)
  # setVariableInferDim!(varid, 0)
  # setVariableInitialized!(vari, false)
  nothing
end

resetVariable!(vari::DFGVariable; solveKey::Symbol=:default  )::Nothing = resetVariable!(getData(vari), solveKey=solveKey)

function resetVariable!(dfg::G,
                        sym::Symbol;
                        solveKey::Symbol=:default  )::Nothing where G <: AbstractDFG
  #
  resetVariable!(getVariable(dfg, sym), solveKey=solveKey)
end




# TODO -- there should be a better way, without retrieving full vertex
# TODO -- Deprecated for DFG -- must update
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

function setDefaultNodeData!(v::DFGVariable,
                             dodims::Int,
                             N::Int,
                             dims::Int;
                             gt=Dict(),
                             initialized::Bool=true,
                             dontmargin::Bool=false,
                             softtype=nothing)::Nothing
  # TODO review and refactor this function, exists as legacy from pre-v0.3.0
  # this should be the only function allocating memory for the node points (unless number of points are changed)
  data = nothing
  if initialized

      pN = AMP.manikde!(randn(dims, N), softtype.manifolds);

    sp = Int[0;] #round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
    gbw = getBW(pN)[:,1]
    gbw2 = Array{Float64}(undef, length(gbw),1)
    gbw2[:,1] = gbw[:]
    pNpts = getPoints(pN)
    #initval, stdev
    setSolverData(v, VariableNodeData(pNpts,
                            gbw2, Symbol[], sp,
                            dims, false, :_null, Symbol[], softtype, true, 0.0, false, dontmargin))
  else
    sp = round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
    setSolverData(v, VariableNodeData(zeros(dims, N),
                            zeros(dims,1), Symbol[], sp,
                            dims, false, :_null, Symbol[], softtype, false, 0.0, false, dontmargin))
  end
  return nothing
end
# if size(initval,2) < N && size(initval, 1) == dims
#   @warn "setDefaultNodeData! -- deprecated use of stdev."
#   p = AMP.manikde!(initval,diag(stdev), softtype.manifolds);
#   pN = resample(p,N)
# if size(initval,2) < N && size(initval, 1) != dims
  # @info "Node value memory allocated but not initialized"
# else
#   pN = AMP.manikde!(initval, softtype.manifolds)
# end
# dims = size(initval,1) # rows indicate dimensions


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
                      smalldata=Dict{String, String}(),
                      checkduplicates::Bool=true  )::DFGVariable where
                        {G <: AbstractDFG,
                         T <: InferenceVariable}
  #
  v = DFGVariable(lbl)
  v.ready = ready
  # v.backendset = backendset
  v.tags = union(labels, Symbol.(softtype.labels), [:VARIABLE])
  v.smallData = smalldata
  setDefaultNodeData!(v, 0, N, softtype.dims, initialized=!autoinit, softtype=softtype, dontmargin=dontmargin) # dodims
  DFG.addVariable!(dfg, v)

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
                      smalldata=Dict{String, String}())::DFGVariable where
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
  @warn "getVal(::Vector{DFGVariable}) is obsolete, use getVal.(DFGVariable) instead."
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
  ftyp = T
  mhcat = parseusermultihypo(multihypo)

  ccw = prepgenericconvolution(Xi, usrfnc, multihypo=mhcat, threadmodel=threadmodel)

  data_ccw = FunctionNodeData{CommonConvWrapper{T}}(Int[], false, false, Int[], Symbol(ftyp.name.module), ccw)
  factor.data = data_ccw

  return factor.data
end

"""
    $SIGNATURES

Returns state of vertex data `.initialized` flag.

Notes:
- used by Bayes tree clique logic.
- similar method in DFG
"""
function isInitialized(vert::Graphs.ExVertex)::Bool
  return getData(vert).initialized
end
# function isInitialized(vert::DFGVariable)::Bool
#   return getData(vert).initialized
# end
# function isInitialized(dfg::T, vsym::Symbol)::Bool where T <: AbstractDFG
#   return isInitialized(DFG.getVariable(dfg, vsym))
# end

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
  varsyms = DFG.getNeighbors(dfg, fct)

  # list of factors to use in init operation
  fctlist = []
  faillist = Symbol[]
  for vsym in varsyms
    xi = DFG.getVariable(dfg, vsym)
    if !isInitialized(xi)
      push!(faillist, vsym)
    end
  end

  # determine if this factor can be used
  canuse = length(varsyms)==1 || (length(faillist)==1 && loovar in faillist)
  if canuse
    push!(fctlist, fct)
  end

  # return if can use, the factor in an array, and the non-initialized variables attached to the factor
  return (canuse, fctlist, faillist )
end


# wow, that was quite far off -- needs testing
# function factorCanInitFromOtherVars(dfg::T,
#                                     fct::Symbol,
#                                     loovar::Symbol)::Tuple{Bool, Vector{Symbol}, Vector{Symbol}} where T <: AbstractDFG
#   #
#   # all variables attached to this factor
#   varsyms = getNeighbors(dfg, fct)
#
#   # list of factors to use in init operation
#   useinitfct = Symbol[]
#   faillist = Symbol[]
#   for vsym in varsyms
#     xi = DFG.getVariable(dfg, vsym)
#     if (isInitialized(xi) && sum(useinitfct .== fct) == 0 ) || length(varsyms) == 1
#       push!(useinitfct, fct)
#     end
#   end
#
#   return (length(useinitfct)==length(varsyms)&&length(faillist)==0,
#           useinitfct,
#           faillist   )
# end

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
                     N::Int=100,
                     logger=ConsoleLogger() )::Bool where T <: AbstractDFG
  #
  didinit = false
  # don't initialize a variable more than once
  if !isInitialized(xi)
    with_logger(logger) do
      @info "try doautoinit! of $(xi.label)"
    end
    # get factors attached to this variable xi
    vsym = xi.label
    neinodes = DFG.getNeighbors(dfg, vsym)
    # proceed if has more than one neighbor OR even if single factor
    if (singles || length(neinodes) > 1)
      # Which of the factors can be used for initialization
      useinitfct = Symbol[]
      # Consider factors connected to $vsym...
      for xifct in neinodes
        canuse, usefct, notusevars = factorCanInitFromOtherVars(dfg, xifct, vsym)
        if canuse
          union!(useinitfct, usefct)
        end
      end
      with_logger(logger) do
        @info "init with useinitfct $useinitfct"
      end
      # println("Consider all singleton (unary) factors to $vsym...")
      # calculate the predicted belief over $vsym
      if length(useinitfct) > 0
        with_logger(logger) do
          @info "do init of $vsym"
        end
        pts,inferdim = predictbelief(dfg, vsym, useinitfct, logger=logger)
        setValKDE!(xi, pts, true, inferdim)
        didinit = true
      end
    end
  end
  return didinit
end

function doautoinit!(dfg::T,
                     Xi::Vector{DFGVariable};
                     singles::Bool=true,
                     N::Int=100,
                     logger=ConsoleLogger() )::Bool where T <: AbstractDFG
  #
  #
  # Mighty inefficient function, since we only need very select fields nearby from a few neighboring nodes
  # do double depth search for variable nodes

  didinit = true

  # loop over all requested variables that must be initialized
  for xi in Xi
    didinit &= doautoinit!(dfg, xi, singles=singles, N=N, logger=logger)
  end
  return didinit
end

function doautoinit!(dfg::T,
                     xsyms::Vector{Symbol};
                     singles::Bool=true,
                     N::Int=100,
                     logger=SimpleLogger(logger)  )::Bool where T <: AbstractDFG
  #
  verts = getVariable.(dfg, xsyms)
  return doautoinit!(dfg, verts, singles=singles, N=N, logger=logger)
end
function doautoinit!(dfg::T,
                     xsym::Symbol;
                     singles::Bool=true,
                     N::Int=100,
                     logger=ConsoleLogger()  )::Bool where T <: AbstractDFG
  #
  return doautoinit!(dfg, [getVariable(dfg, xsym);], singles=singles, N=N, logger=logger)
end

"""
    $(SIGNATURES)

Workaround function when first-version (factor graph based) auto initialization fails.  Usually occurs when using factors that have high connectivity to multiple variables.
"""
function manualinit!(dfg::T, vert::DFGVariable, pX::BallTreeDensity)::Nothing where T <: AbstractDFG
  setValKDE!(vert, pX, true)
  # getData(vert).initialized = true
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
  setValKDE!(vert, Xpre, true) # dfg, sym
  # getData(dfg, sym).initialized = true
  return nothing
end
function manualinit!(dfg::G, sym::Symbol, pts::Array{Float64,2}) where G <: AbstractDFG
  var = getVariable(dfg, sym)
  pp = manikde!(pts, getManifolds(var))
  manualinit!(dfg,sym,pp)
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

function assembleFactorName(dfg::T, Xi::Vector{DFGVariable}; maxparallel::Int=50)::Symbol where T <: AbstractDFG
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
  return Symbol(namestring)
end

"""
    $(SIGNATURES)

Add factor with user defined type <: FunctorInferenceType to the factor graph
object. Define whether the automatic initialization of variables should be
performed.  Use order sensitive `multihypo` keyword argument to define if any
variables are related to data association uncertainty.
"""
function addFactor!(dfg::G,
                    Xi::Vector{DFGVariable},
                    usrfnc::R;
                    multihypo::Union{Nothing,Tuple,Vector{Float64}}=nothing,
                    ready::Int=1,
                    labels::Vector{Symbol}=Symbol[],
                    autoinit::Bool=true,
                    threadmodel=SingleThreaded,
                    maxparallel::Int=50  ) where
                      {G <: AbstractDFG,
                       R <: Union{FunctorInferenceType, InferenceType}}
  #
  namestring = assembleFactorName(dfg, Xi, maxparallel=maxparallel)
  newFactor = DFGFactor{CommonConvWrapper{R}, Symbol}(Symbol(namestring))
  newFactor.tags = union(labels, [:FACTOR]) # TODO: And session info
  # addNewFncVertInGraph!(fgl, newvert, currid, namestring, ready)
  newData = setDefaultFactorNode!(dfg, newFactor, Xi, deepcopy(usrfnc), multihypo=multihypo, threadmodel=threadmodel)

  # TODO: Need to remove this...
  for vert in Xi
    push!(newData.fncargvID, vert.label) # vert._internalId # YUCK :/ -- Yup, this is a problem
  end

  success = DFG.addFactor!(dfg, Xi, newFactor)

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
      threadmodel=SingleThreaded,
      maxparallel::Int=50  ) where
        {G <: AbstractDFG,
         R <: Union{FunctorInferenceType, InferenceType}}
  #
  verts = map(vid -> DFG.getVariable(dfg, vid), xisyms)
  addFactor!(dfg, verts, usrfnc, multihypo=multihypo, ready=ready, labels=labels, autoinit=autoinit, threadmodel=threadmodel, maxparallel=maxparallel )
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

"""
    $SIGNATURES

Determine the variable ordering used to construct both the Bayes Net and Bayes/Junction/Elimination tree.

Notes
- Heuristic method -- equivalent to QR or Cholesky.
- Are using Blas `QR` function to extract variable ordering.
- **NOT USING SUITE SPARSE** -- which would requires commercial license.
- For now `A::Array{<:Number,2}` as a dense matrix.
- Columns of `A` are system variables, rows are factors (without differentiating between partial or full factor).

Future
- TODO: `A` should be sparse data structure (when we exceed 10'000 var dims)
"""
function getEliminationOrder(dfg::G; ordering::Symbol=:qr) where G <: AbstractDFG
  # Get the sparse adjacency matrix, variable, and factor labels
  adjMat, permuteds, permutedsf = DFG.getAdjacencyMatrixSparse(dfg)

  # Create dense adjacency matrix
  A = Array(adjMat)

  p = Int[]
  if ordering==:chol
      p = cholfact(A'A,:U,Val(true))[:p] #,pivot=true
  elseif ordering==:qr
    # this is the default
    q,r,p = qr(A, Val(true))
  else
    prtslperr("getEliminationOrder -- cannot do the requested ordering $(ordering)")
  end

  # Return the variable ordering that we should use for the Bayes map
  return permuteds[p]
end


# lets create all the vertices first and then deal with the elimination variables thereafter
function addBayesNetVerts!(dfg::G, elimOrder::Array{Symbol,1}) where G <: AbstractDFG
  for pId in elimOrder
    vert = DFG.getVariable(dfg, pId)
    if getData(vert).BayesNetVertID == nothing || getData(vert).BayesNetVertID == :_null # Special serialization case of nothing
      @debug "[AddBayesNetVerts] Assigning $pId.data.BayesNetVertID = $pId"
      getData(vert).BayesNetVertID = pId
    else
      @warn "addBayesNetVerts -- Something is wrong, variable '$pId' should not have an existing Bayes net reference to '$(getData(vert).BayesNetVertID)'"
    end
  end
end

function addConditional!(dfg::G, vertId::Symbol, Si::Vector{Symbol})::Nothing where G <: AbstractDFG
  bnv = DFG.getVariable(dfg, vertId)
  bnvd = getData(bnv) # bnv.attributes["data"]
  bnvd.separator = Si
  for s in Si
    push!(bnvd.BayesNetOutVertIDs, s)
  end
  return nothing
end

function addChainRuleMarginal!(dfg::G, Si::Vector{Symbol}; maxparallel::Int=50) where G <: AbstractDFG
    @show Si
  lbls = String[]
  genmarg = GenericMarginal()
  Xi = map(v -> DFG.getVariable(dfg, v), Si)
  # @info "adding marginal to"
  # for x in Xi
  #   @info "x.index=",x.index
  # end
  addFactor!(dfg, Xi, genmarg, autoinit=false, maxparallel=maxparallel)
  nothing
end

function rmVarFromMarg(dfg::G,
                       fromvert::DFGVariable,
                       gm::Vector{DFGFactor};
                       maxparallel::Int=50 )::Nothing where G <: AbstractDFG
  @info " - Removing $(fromvert.label)"
  for m in gm
    @info "Looking at $(m.label)"
    for n in DFG.getNeighbors(dfg, m) #x1, x2
      if n == fromvert.label # n.label ==? x1
        @info "   - Breaking link $(m.label)->$(fromvert.label)..."
        @info "     - Original links: $(DFG.ls(dfg, m))"
        remvars = setdiff(DFG.ls(dfg, m), [fromvert.label])
        @info "     - New links: $remvars"

        DFG.deleteFactor!(dfg, m) # Remove it
        if length(remvars) > 0
          @info "$(m.label) still has links to other variables, readding it back..."
          addFactor!(dfg, remvars, getData(m).fnc.usrfnc!, autoinit=false, maxparallel=maxparallel)
        else
          @info "$(m.label) doesn't have any other links, not adding it back..."
        end
      end
    end
    # Added back in chain rule.
    if DFG.exists(dfg, m) && length(DFG.getNeighbors(dfg, m)) <= 1
      @warn "removing vertex id=$(m.label)"
      DFG.deleteFactor!(dfg, m)
    end
  end
  return nothing
end

function buildBayesNet!(dfg::G,
                        elimorder::Vector{Symbol};
                        maxparallel::Int=50)::Nothing where G <: AbstractDFG
  #
  # addBayesNetVerts!(dfg, elimorder)
  for v in elimorder
    @info ""
    @info "Eliminating $(v)"
    @info "==============="
    @info ""
    # which variable are we eliminating

    # all factors adjacent to this variable
    fi = Symbol[]
    Si = Symbol[]
    gm = DFGFactor[]

    vert = DFG.getVariable(dfg, v)
    for fctId in DFG.getNeighbors(dfg, vert)
      fct = DFG.getFactor(dfg, fctId)
      if (getData(fct).eliminated != true)
        push!(fi, fctId)
        for sepNode in DFG.getNeighbors(dfg, fct)
          # TODO -- validate !(sepNode.index in Si) vs. older !(sepNode in Si)
          if sepNode != v && !(sepNode in Si) # Symbol comparison!
            push!(Si,sepNode)
          end
        end
        getData(fct).eliminated = true #fct.attributes["data"].eliminated = true
      end

      if typeof(getData(fct).fnc) == CommonConvWrapper{GenericMarginal}
        push!(gm, fct)
      end
    end

    if v != elimorder[end]
      addConditional!(dfg, v, Si)
      # not yet inserting the new prior p(Si) back into the factor graph
    end

    # mark variable
    getData(vert).eliminated = true

    # TODO -- remove links from current vertex to any marginals
    rmVarFromMarg(dfg, vert, gm, maxparallel=maxparallel)

    #add marginal on remaining variables... ? f(xyz) = f(x | yz) f(yz)
    # new function between all Si (round the outside, right the outside)
    length(Si) > 0 && addChainRuleMarginal!(dfg, Si, maxparallel=maxparallel)

  end
  return nothing
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

function getKDE(vnd::VariableNodeData)
  AMP.manikde!(getVal(vnd), getBW(vnd)[:,1], getSofttype(vnd).manifolds)
end


"""
    $(SIGNATURES)

Get KernelDensityEstimate kde estimate stored in variable node.
"""
function getKDE(v::DFGVariable)
  return getKDE(getData(v))
end

function getVert(dfg::G, sym::Symbol, nt::Symbol=:var) where G <: AbstractDFG
  @warn "IIF.getVert is deprecated, use DFG.getVariable or DFG.getFactor instead."
  if nt == :var
    return DFG.getVariable(dfg, sym)
  elseif nt == :fct
    return DFG.getFactor(dfg, sym)
  else
    error("unknown getVert request nt=$nt")
  end
end

"""
    $(SIGNATURES)

Get KernelDensityEstimate kde estimate stored in variable node.
"""
function getVertKDE(v::DFGVariable)
  return getKDE(v)
end
function getVertKDE(dfg::G, id::Int) where G <: AbstractDFG
  v = DFG.getVariable(dfg, id)
  return getKDE(v)
end
function getVertKDE(dfg::G, lbl::Symbol) where G <: AbstractDFG
  v = DFG.getVariable(dfg, lbl)
  return getKDE(v)
end
function getKDE(dfg::G, lbl::Symbol) where G <: AbstractDFG
  return getVertKDE(dfg, lbl)
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
