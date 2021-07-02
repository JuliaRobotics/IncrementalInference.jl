


"""
$SIGNATURES

Initialize an empty in-memory DistributedFactorGraph `::DistributedFactorGraph` object.
"""
function initfg(dfg::T=InMemDFGType(solverParams=SolverParams());
                                sessionname="NA",
                                robotname="",
                                username="",
                                cloudgraph=nothing) where T <: AbstractDFG
#
  #
  return dfg
end


#init an empty fg with a provided type and SolverParams
function initfg(::Type{T};solverParams=SolverParams(),
                      sessionname="NA",
                      robotname="",
                      username="",
                      cloudgraph=nothing) where T <: AbstractDFG
  #
  return T(solverParams=solverParams)
end

function initfg(::Type{T},solverParams::SolverParams;
                      sessionname="NA",
                      robotname="",
                      username="",
                      cloudgraph=nothing) where T <: AbstractDFG
  #
  return T{SolverParams}(solverParams=solverParams)
end



# Should deprecate in favor of TensorCast.jl
reshapeVec2Mat(vec::Vector, rows::Int) = reshape(vec, rows, round(Int,length(vec)/rows))



## ==============================================================================================
## MOVE TO / CONSOLIDATE WITH DFG
## ==============================================================================================


"""
    $(SIGNATURES)

Fetch the variable marginal joint sampled points.  Use [`getBelief`](@ref) to retrieve the full Belief object.
"""
getVal(v::DFGVariable; solveKey::Symbol=:default) = v.solverDataDict[solveKey].val
getVal(v::DFGVariable, idx::Int; solveKey::Symbol=:default) = v.solverDataDict[solveKey].val[:,idx]
getVal(vnd::VariableNodeData) = vnd.val
getVal(vnd::VariableNodeData, idx::Int) = vnd.val[:, idx]
getVal(dfg::AbstractDFG, lbl::Symbol; solveKey::Symbol=:default) = getVariable(dfg, lbl).solverDataDict[solveKey].val

"""
    $(SIGNATURES)

Get the number of points used for the current marginal belief estimate represtation for a particular variable in the factor graph.
"""
function getNumPts(v::DFGVariable; solveKey::Symbol=:default)::Int
  return length(getVal(getSolverData(v, solveKey)))
end

function getBW(vnd::VariableNodeData)
  return vnd.bw
end

# setVal! assumes you will update values to database separate, this used for local graph mods only
function getBWVal(v::DFGVariable; solveKey::Symbol=:default)
  return getSolverData(v, solveKey).bw
end
function setBW!(vd::VariableNodeData, bw::Array{Float64,2}; solveKey::Symbol=:default)
  vd.bw = bw
  nothing
end
function setBW!(v::DFGVariable, bw::Array{Float64,2}; solveKey::Symbol=:default)
  setBW!(getSolverData(v, solveKey), bw)
  nothing
end

function setVal!(vd::VariableNodeData, val::AbstractVector{P}) where P
  vd.val = val
  nothing
end
function setVal!(v::DFGVariable, val::AbstractVector{P}; solveKey::Symbol=:default) where P
    setVal!(getSolverData(v, solveKey), val)
    nothing
end
function setVal!(vd::VariableNodeData, val::AbstractVector{P}, bw::Array{Float64,2}) where P
    setVal!(vd, val)
    setBW!(vd, bw)
    nothing
end
function setVal!(v::DFGVariable, val::AbstractVector{P}, bw::Array{Float64,2}; solveKey::Symbol=:default) where P
  setVal!(v, val, solveKey=solveKey)
  setBW!(v, bw, solveKey=solveKey)
  nothing
end
function setVal!(vd::VariableNodeData, val::AbstractVector{P}, bw::Vector{Float64}) where P
  setVal!(vd, val, reshape(bw,length(bw),1))
  nothing
end
function setVal!(v::DFGVariable, val::AbstractVector{P}, bw::Vector{Float64}; solveKey::Symbol=:default) where P
  setVal!(getSolverData(v, solveKey=solveKey), val, bw)
  nothing
end
function setVal!(dfg::AbstractDFG, sym::Symbol, val::AbstractVector{P}; solveKey::Symbol=:default) where P
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
                    pts::AbstractVector{P},
                    bws::Vector{Float64},
                    setinit::Bool=true,
                    inferdim::Float64=0.0 ) where P
  #
  setVal!(vd, pts, bws) # BUG ...al!(., val, . ) ## TODO -- this can be a little faster
  setinit ? (vd.initialized = true) : nothing
  vd.inferdim = inferdim
  nothing
end

function setValKDE!(vd::VariableNodeData,
                    val::AbstractVector{P},
                    setinit::Bool=true,
                    inferdim::Real=0.0  ) where P
  # recover variableType information
  varType = getVariableType(vd)
  p = AMP.manikde!(val, varType)
  setValKDE!(vd, p, setinit, inferdim)
  nothing
end

function setValKDE!(v::DFGVariable,
                    val::AbstractVector{P},
                    bws::Array{<:Real,2},
                    setinit::Bool=true,
                    inferdim::Float64=0;
                    solveKey::Symbol=:default) where P
  # recover variableType information
  setValKDE!(getSolverData(v, solveKey), val, bws[:,1], setinit, inferdim )

  nothing
end

function setValKDE!(v::DFGVariable,
                    val::AbstractVector{P},
                    setinit::Bool=true,
                    inferdim::Float64=0.0;
                    solveKey::Symbol=:default) where P
  # recover variableType information
  setValKDE!(getSolverData(v, solveKey),val, setinit, inferdim )
  nothing
end
function setValKDE!(v::DFGVariable,
                    em::TreeBelief,
                    setinit::Bool=true;
                    # inferdim::Union{Float32, Float64, Int32, Int64}=0;
                    solveKey::Symbol=:default  )
  #
  setValKDE!(v, em.val, em.bw, setinit, em.inferdim, solveKey=solveKey)
  nothing
end
function setValKDE!(v::DFGVariable,
                    p::ManifoldKernelDensity,
                    setinit::Bool=true,
                    inferdim::Union{Float32, Float64, Int32, Int64}=0;
                    solveKey::Symbol=:default  )
  #
  # @error("TESTING setValKDE! ", solveKey, string(listSolveKeys(v)))
  setValKDE!(getSolverData(v,solveKey),p,setinit,Float64(inferdim))
  nothing
end
function setValKDE!(dfg::AbstractDFG,
                    sym::Symbol,
                    p::ManifoldKernelDensity,
                    setinit::Bool=true,
                    inferdim::Union{Float32, Float64, Int32, Int64}=0;
                    solveKey::Symbol=:default  )
  #
  setValKDE!(getVariable(dfg, sym), p, setinit, inferdim, solveKey=solveKey)
  nothing
end


"""
    $SIGNATURES

Set variable initialized status.
"""
function setVariableInitialized!( varid::VariableNodeData,
                                  status::Bool)
  #
  varid.initialized = status
end

setVariableInitialized!(vari::DFGVariable, status::Bool) = setVariableInitialized!(getSolverData(vari), status)


"""
    $SIGNATURES

Set method for the inferred dimension value in a variable.
"""
setVariableInferDim!(varid::VariableNodeData, val::Real) = varid.inferdim = convert(Float64,val)
setVariableInferDim!(vari::DFGVariable, val::Real) = setVariableInferDim!(getSolverData(vari), val)

## ==============================================================================================
## ==============================================================================================


function setValKDE!(vd::VariableNodeData,
                    p::ManifoldKernelDensity,
                    setinit::Bool=true,
                    inferdim::Union{Float32, Float64, Int32, Int64}=0 )
  #
  ptsArr = AMP.getPoints(p)
  # @show typeof(ptsArr)
  # @cast ptsArr[j][i] := pts[i,j]
  bws = getBW(p)[:,1]
  setValKDE!(vd,ptsArr,bws,setinit,inferdim )
  nothing
end

"""
    $SIGNATURES

Reset the solve state of a variable to uninitialized/unsolved state.
"""
function resetVariable!(varid::VariableNodeData;
                        solveKey::Symbol=:default  )::Nothing
  #
  val = getBelief(varid)
  pts = AMP.getPoints(val)
  # TODO not all manifolds will initialize to zero
  for pt in pts
    fill!(pt, 0.0)
  end
  pn = manikde!(pts, zeros(AMP.Ndim(val)), getManifolds(varid))
  setValKDE!(varid, pn, false, 0.0)
  # setVariableInferDim!(varid, 0)
  # setVariableInitialized!(vari, false)
  nothing
end

resetVariable!(vari::DFGVariable; solveKey::Symbol=:default  )::Nothing = resetVariable!(getSolverData(vari), solveKey=solveKey)

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



function DefaultNodeDataParametric( dodims::Int,
                                    dims::Int,
                                    variableType::InferenceVariable;
                                    initialized::Bool=true,
                                    dontmargin::Bool=false)::VariableNodeData

  # this should be the only function allocating memory for the node points
  if false && initialized
    error("not implemented yet")
    # pN = AMP.manikde!(randn(dims, N), variableType.manifolds);
    #
    # sp = Int[0;] #round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
    # gbw = getBW(pN)[:,1]
    # gbw2 = Array{Float64}(undef, length(gbw),1)
    # gbw2[:,1] = gbw[:]
    # pNpts = getPoints(pN)
    # #initval, stdev
    # return VariableNodeData(pNpts,
    #                         gbw2, Symbol[], sp,
    #                         dims, false, :_null, Symbol[], variableType, true, 0.0, false, dontmargin)
  else
    sp = round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
    return VariableNodeData([zeros(dims) for _ in 1:1],
                            zeros(dims,dims), Symbol[], sp,
                            dims, false, :_null, Symbol[], variableType, false, 0.0, false, dontmargin, 0, 0, :parametric)
  end

end

function setDefaultNodeDataParametric!(v::DFGVariable, variableType::InferenceVariable; kwargs...)
  vnd = DefaultNodeDataParametric(0, variableType |> getDimension, variableType; kwargs...)
  setSolverData!(v, vnd, :parametric)
  return nothing
end

"""
    $SIGNATURES

Create new solverData.

Notes
- Used during creation of new variable, as well as in CSM unique `solveKey`.
"""
function setDefaultNodeData!( v::DFGVariable,
                              dodims::Int,
                              N::Int,
                              dims::Int;
                              solveKey::Symbol=:default,
                              gt=Dict(),
                              initialized::Bool=true,
                              dontmargin::Bool=false,
                              varType=nothing )
  #
  # TODO review and refactor this function, exists as legacy from pre-v0.3.0
  # this should be the only function allocating memory for the node points (unless number of points are changed)
  data = nothing
  if initialized

      pN = AMP.manikde!(randn(dims, N), getManifolds(varType));

    sp = Int[0;] #round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
    gbw = getBW(pN)[:,1]
    gbw2 = Array{Float64}(undef, length(gbw),1)
    gbw2[:,1] = gbw[:]
    pNpts = getPoints(pN)
    #initval, stdev
    setSolverData!(v, VariableNodeData(pNpts,
                            gbw2, Symbol[], sp,
                            dims, false, :_null, Symbol[], 
                            varType, true, 0.0, false, dontmargin,0,0,solveKey), solveKey)
  else
    sp = round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
    valpts = Vector{getPointType(varType)}(undef,N)
    for i in 1:length(valpts)
      valpts[i] = getPointIdentity(varType)
    end
    bws = zeros(dims,1)
    setSolverData!(v, VariableNodeData(valpts, bws,
                                        Symbol[], sp,
                                        dims, false, :_null, Symbol[], 
                                        varType, false, 0.0, false, dontmargin,0,0,solveKey), solveKey)
    #
  end
  return nothing
end
# if size(initval,2) < N && size(initval, 1) == dims
#   @warn "setDefaultNodeData! -- deprecated use of stdev."
#   p = AMP.manikde!(initval,diag(stdev), varType.manifolds);
#   pN = resample(p,N)
# if size(initval,2) < N && size(initval, 1) != dims
  # @info "Node value memory allocated but not initialized"
# else
#   pN = AMP.manikde!(initval, varType.manifolds)
# end
# dims = size(initval,1) # rows indicate dimensions



"""
    $SIGNATURES

Reference data can be stored in the factor graph as a super-solve.

Notes
- Intended as a mechanism to store reference data alongside the numerical computations.
"""
function setVariableRefence!( dfg::AbstractDFG,
                              sym::Symbol,
                              val::AbstractVector;
                              refKey::Symbol=:reference)
  #
  # which variable to update
  var = getVariable(dfg, sym)

  # Construct an empty VND object
  vnd = VariableNodeData( val,
                          zeros(getDimension(var),1),
                          Symbol[],
                          Int[0;],
                          getDimension(var),
                          false,
                          :_null,
                          Symbol[],
                          getVariableType(var),
                          true,
                          0.0,
                          false,
                          true  )
  #
  # set the value in the DFGVariable
  setSolverData!(var, vnd, refKey)
end


# get instance from variableType
_variableType(varType::InferenceVariable) = varType
_variableType(varType::Type{<:InferenceVariable}) = varType()

"""
$(SIGNATURES)

Add a variable node `label::Symbol` to `dfg::AbstractDFG`, as `varType<:InferenceVariable`.

Notes
-----
- keyword `nanosecondtime` is experimental and intended as the whole subsection portion -- i.e. accurateTime = (timestamp MOD second) + Nanosecond

Example
-------

```julia
fg = initfg()
addVariable!(fg, :x0, Pose2)
```
"""
function addVariable!(dfg::AbstractDFG,
                      label::Symbol,
                      varTypeU::Union{T, Type{T}}; 
                      N::Int=getSolverParams(dfg).N,
                      solvable::Int=1,
                      timestamp::Union{DateTime,ZonedDateTime}=now(localzone()),
                      nanosecondtime::Union{Nanosecond,Int64,Nothing}=nothing,
                      dontmargin::Bool=false,
                      labels::Union{Vector{Symbol},Nothing}=nothing,
                      tags::Vector{Symbol}=Symbol[],
                      smalldata=Dict{Symbol, DFG.SmallDataTypes}(),
                      checkduplicates::Bool=true,
                      initsolvekeys::Vector{Symbol}=getSolverParams(dfg).algorithms ) where T<:InferenceVariable
  #
  varType = _variableType(varTypeU)
  # TODO Remove deprecation in v0.16 
  if :ut in fieldnames(T)	
    Base.depwarn("Field `ut` (microseconds) for variable type ($T) has been deprecated please use DFGVariable.nstime, kwarg: nanosecondtime", :addVariable!)	
    if isnothing(nanosecondtime)
      varType.ut == -9999999999 && error("please define a time for type $(T), use FGVariable.nstime, kwarg: nanosecondtime")
      nanosecondtime = Nanosecond(varType.ut*1000)	
    else 	
      @warn "Nanosecond time has been specified as $nanosecondtime, ignoring `ut` field value: $(varType.ut)."	
    end	
  elseif isnothing(nanosecondtime)	
    nanosecondtime = Nanosecond(0)	
  end

  # deprecate in v0.16
  labels isa Vector ? (union!(tags, labels); @warn("labels is deprecated, use tags instead")) : nothing
  
  union!(tags, [:VARIABLE])
  v = DFGVariable(label, varType; tags=Set(tags), smallData=smalldata, solvable=solvable, timestamp=timestamp, nstime=Nanosecond(nanosecondtime))

  (:default in initsolvekeys) &&
    setDefaultNodeData!(v, 0, N, getDimension(varType), initialized=false, varType=varType, dontmargin=dontmargin) # dodims

  (:parametric in initsolvekeys) &&
    setDefaultNodeDataParametric!(v, varType, initialized=false, dontmargin=dontmargin)

  DFG.addVariable!(dfg, v)

  return v
end

function _resizePointsVector!(vecP::AbstractVector{P}, mkd::ManifoldKernelDensity, N::Int) where P
  #
  pN = length(vecP)
  resize!(vecP, N)
  for j in pN:N
    smp = AMP.sample(mkd, 1)[1]
    @show j, smp, typeof(smp), typeof(vecP[j])
    vecP[j] = smp
  end

  vecP
end


"""
    $(SIGNATURES)

Prepare the particle arrays `ARR` to be used for approximate convolution.
This function ensures that ARR has te same dimensions among all the parameters.
Function returns with ARR[sfidx] pointing at newly allocated deepcopy of the
existing values in getVal(Xi[.label==solvefor]).

Notes
- Return values `sfidx` is the element in ARR where `Xi.label==solvefor` and
- `maxlen` is length of all (possibly resampled) `ARR` contained particles.
- `Xi` is order sensitive.
- for initialization, solveFor = Nothing.
- `P = getPointType(<:InferenceVariable)`
"""
function prepareparamsarray!( ARR::AbstractVector{<:AbstractVector{P}},
                              Xi::Vector{<:DFGVariable},
                              solvefor::Union{Nothing, Symbol},
                              N::Int=0;
                              solveKey::Symbol=:default  ) where P
  #
  LEN = Int[]
  maxlen = N # FIXME see #105
  count = 0
  sfidx = 0

  for xi in Xi
    vecP = getVal(xi, solveKey=solveKey)
    push!(ARR, vecP)
    LEN = length.(ARR)
    maxlen = maximum(LEN)
    count += 1
    if xi.label == solvefor
      sfidx = count #xi.index
    end
  end

  # resample variables with too few kernels (manifolds points)
  SAMP = LEN .< maxlen
  for i in 1:count
    if SAMP[i]
      Pr = getBelief(Xi[i], solveKey)
      _resizePointsVector!(ARR[i], Pr, maxlen)
    end
  end

  # TODO --rather define reusable memory for the proposal
  # we are generating a proposal distribution, not direct replacement for existing memory and hence the deepcopy.
  if sfidx > 0 
    ARR[sfidx] = deepcopy(ARR[sfidx]) 
  end

  # get solvefor manifolds
  mani = length(Xi)==0 || sfidx==0 ? (:null,) : getManifold(Xi[sfidx])

  # FIXME, forcing maxlen to N results in errors (see test/testVariousNSolveSize.jl) see #105
  # maxlen = N == 0 ? maxlen : N
  return maxlen, sfidx, mani
end

function parseusermultihypo(multihypo::Nothing, nullhypo::Float64)
  verts = Symbol[]
  mh = nothing
  return mh, nullhypo
end
function parseusermultihypo(multihypo::Vector{Float64}, nullhypo::Float64)
  mh = nothing
  if 0 < length(multihypo)
    multihypo2 = multihypo
    multihypo2[1-1e-10 .< multihypo] .= 0.0
    # check that terms sum to full probability
    @assert abs(sum(multihypo2) % 1) < 1e-10  || 1-1e-10 < sum(multihypo2) % 1 "ensure multihypo sums to a (or nearly, 1e-10) interger, see #1086"
    # check that only one variable broken into fractions
    @assert sum(multihypo2[1e-10 .< multihypo2]) ≈ 1
    # force normalize something that is now known to be close
    multihypo2 ./= sum(multihypo2)

    mh = Categorical(Float64[multihypo2...] )
  end
  return mh, nullhypo
end

# import IncrementalInference: prepgenericconvolution, convert

"""
    $SIGNATURES

Function to calculate measurement dimension from factor sampling.

Notes
- Will not work in all situations, but good enough so far.
  - # TODO standardize via domain or manifold definition...??
"""
function calcZDim(cf::CalcFactor{T}) where {T <: FunctorInferenceType}
  #
  # zdim = T != GenericMarginal ? size(getSample(usrfnc, 2)[1],1) : 0
  zdim = if T != GenericMarginal
    # vnds = Xi # (x->getSolverData(x)).(Xi)
    # NOTE try to make sure we get matrix back (not a vector)
    smpls = sampleFactor(cf, 2)[1]
    length(smpls[1])
  else
    0
  end
  return zdim
end


function prepgenericconvolution(Xi::Vector{<:DFGVariable},
                                usrfnc::T;
                                multihypo::Union{Nothing, Distributions.Categorical}=nothing,
                                nullhypo::Real=0.0,
                                threadmodel=MultiThreaded,
                                inflation::Real=0.0  ) where {T <: FunctorInferenceType}
  #
  pttypes = getVariableType.(Xi) .|> getPointType
  PointType = 0 < length(pttypes) ? pttypes[1] : Vector{Float64}
  ARR = Vector{Vector{PointType}}()
  maxlen, sfidx, mani = prepareparamsarray!(ARR, Xi, nothing, 0) # Nothing for init.

  # standard factor metadata
  sflbl = 0==length(Xi) ? :null : getLabel(Xi[end])
  fmd = FactorMetadata(Xi, getLabel.(Xi), ARR, sflbl, nothing)
  # guess measurement points type
  MeasType = Vector{Float64} # FIXME use `usrfnc` to get this information instead
  cf = CalcFactor( usrfnc, fmd, 0, 1, (Vector{MeasType}(),), ARR)

  zdim = calcZDim(cf)
  # zdim = T != GenericMarginal ? size(getSample(usrfnc, 2)[1],1) : 0
  certainhypo = multihypo !== nothing ? collect(1:length(multihypo.p))[multihypo.p .== 0.0] : collect(1:length(Xi))
  
  # sort out partialDims here
  ispartl = hasfield(T, :partial)
  partialDims = if ispartl
    Int[usrfnc.partial...]
  else
    Int[]
  end

  ccw = CommonConvWrapper(
          usrfnc,
          PointType[],
          zdim,
          ARR,
          fmd,
          specialzDim = hasfield(T, :zDim),
          partial = ispartl,
          hypotheses=multihypo,
          certainhypo=certainhypo,
          nullhypo=nullhypo,
          threadmodel=threadmodel,
          inflation=inflation,
          partialDims=partialDims
        )
  #
  return ccw
end

# TODO perhaps consolidate with constructor?
"""
$SIGNATURES

Generate the default factor data for a new DFGFactor.
"""
function getDefaultFactorData(dfg::AbstractDFG,
                              Xi::Vector{<:DFGVariable},
                              usrfnc::T;
                              multihypo::Vector{<:Real}=Float64[],
                              nullhypo::Float64=0.0,
                              threadmodel=SingleThreaded,
                              eliminated::Bool = false,
                              potentialused::Bool = false,
                              edgeIDs = Int[],
                              solveInProgress = 0,
                              inflation::Real=getSolverParams(dfg).inflation ) where T <: FunctorInferenceType
  #

  # prepare multihypo particulars
  # storeMH::Vector{Float64} = multihypo == nothing ? Float64[] : [multihypo...]
  mhcat, nh = parseusermultihypo(multihypo, nullhypo)

  # allocate temporary state for convolutional operations (not stored)
  ccw = prepgenericconvolution(Xi, usrfnc, multihypo=mhcat, nullhypo=nh, threadmodel=threadmodel, inflation=inflation)

  # and the factor data itself
  return FunctionNodeData{CommonConvWrapper{T}}(eliminated, potentialused, edgeIDs, ccw, multihypo, ccw.certainhypo, nullhypo, solveInProgress, inflation)
end


"""
    $SIGNATURES

Return `::Bool` on whether at least one hypothesis is available for intended computations (assuming direction `sfidx`).
"""
function isLeastOneHypoAvailable( sfidx::Int,
                                  certainidx::Vector{Int},
                                  uncertnidx::Vector{Int},
                                  isinit::Vector{Bool})
  #
  # @show isinit
  # @show sfidx in certainidx, sum(isinit[uncertnidx])
  # @show sfidx in uncertnidx, sum(isinit[certainidx])
  return  sfidx in certainidx && 0 < sum(isinit[uncertnidx]) ||
          sfidx in uncertnidx && sum(isinit[certainidx]) == length(certainidx)
end

"""
    $SIGNATURES

Return `(::Bool, ::OKVarlist, ::NotOkayVarList)` on whether all other variables (besides `loovar::Symbol`)
attached to factor `fct::Symbol` are all initialized -- i.e. `fct` is usable.

Notes:
- Special carve out for multihypo cases, see issue 427, where at least one hypothesis should be available, but not all required at first.

Development Notes
* TODO get faster version of isInitialized for database version

Related

doautoinit!, initManual!, isInitialized, isMultihypo
"""
function factorCanInitFromOtherVars(dfg::AbstractDFG,
                                    fct::Symbol,
                                    loovar::Symbol;
                                    solveKey::Symbol=:default)
  #
  # all variables attached to this factor
  varsyms = DFG.getNeighbors(dfg, fct)

  # which element is being solved for
  sfidx = (1:length(varsyms))[varsyms .== loovar][1]
  # list of factors to use in init operation
  fctlist = Symbol[]
  # list fo variables that cannot be used
  faillist = Symbol[]
  isinit = Bool[]
  for vsym in varsyms
    # check each variable one by one
    xi = DFG.getVariable(dfg, vsym)
    isi = isInitialized(xi, solveKey)
    push!(isinit, isi)
    if !isi
      push!(faillist, vsym)
    end
  end

  ## determine if this factor can be used
  # priors and general n-ary cases
  canuse = length(varsyms)==1 || (length(faillist)==1 && loovar in faillist)
  ## special multihypo case (at least one hypothesis is available or initializing first hypo)
  fctnode = getFactor(dfg, fct)
  # @show canuse, isMultihypo(fctnode), isinit
  if !canuse && isMultihypo(fctnode)
    # multihypo=[1;0.5;0.5] : sfidx=1, isinit=[0,1,0] -- true
    # multihypo=[1;0.5;0.5] : sfidx=1, isinit=[0,0,1] -- true
    # multihypo=[1;0.5;0.5] : sfidx=2|3, isinit=[1,0,0] -- true
    mhp = getMultihypoDistribution(fctnode).p
    allmhp,certainidx,uncertnidx = getHypothesesVectors(mhp)
    if isLeastOneHypoAvailable(sfidx, certainidx, uncertnidx, isinit)
       # special case works
       @info "allowing init from incomplete set of previously initialized hypotheses, fct=$fct"
       canuse = true
    end
  end

  # should add the factor for use?
  if canuse
    push!(fctlist, fct)
  end

  # return if can use, the factor in an array, and the non-initialized variables attached to the factor
  return (canuse, fctlist, faillist)::Tuple{Bool, Vector{Symbol}, Vector{Symbol}}
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
factor graph `fgl`.  Possibly called from `addFactor!`, or `doCliqAutoInitUp!` (?).

Notes:
- Special carve out for multihypo cases, see issue 427.

Development Notes:
> Target factor is first (singletons) or second (dim 2 pairwise) variable vertex in `xi`.
* TODO use DFG properly with local operations and DB update at end.
* TODO get faster version of `isInitialized` for database version.
* TODO: Persist this back if we want to here.
"""
function doautoinit!( dfg::AbstractDFG,
                      xi::DFGVariable;
                      solveKey::Symbol=:default,
                      singles::Bool=true,
                      N::Int=100,
                      logger=ConsoleLogger() )
  #
  didinit = false
  # don't initialize a variable more than once
  if !isInitialized(xi, solveKey)
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
        canuse, usefct, notusevars = factorCanInitFromOtherVars(dfg, xifct, vsym, solveKey=solveKey)
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
        # FIXME ensure a product of only partial densities and returned pts are put to proper dimensions
        pts,inferdim = predictbelief(dfg, vsym, useinitfct, solveKey=solveKey, logger=logger)
        setValKDE!(xi, pts, true, inferdim, solveKey=solveKey)
        # Update the estimates (longer DFG function used so cloud is also updated)
        setVariablePosteriorEstimates!(dfg, xi.label, solveKey)
        # Update the data in the event that it's not local
        # TODO perhaps usecopy=false
        updateVariableSolverData!(dfg, xi, solveKey, true; warn_if_absent=false)    
        # deepcopy graphinit value, see IIF #612
        updateVariableSolverData!(dfg, xi.label, getSolverData(xi, solveKey), :graphinit, true, Symbol[]; warn_if_absent=false)
        didinit = true
      end
    end
  end
  return didinit
end

function doautoinit!( dfg::T,
                      Xi::Vector{<:DFGVariable};
                      solveKey::Symbol=:default,
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
    didinit &= doautoinit!(dfg, xi, solveKey=solveKey, singles=singles, N=N, logger=logger)
  end
  return didinit
end

function doautoinit!( dfg::T,
                      xsyms::Vector{Symbol};
                      solveKey::Symbol=:default,
                      singles::Bool=true,
                      N::Int=100,
                      logger=ConsoleLogger()  )::Bool where T <: AbstractDFG
  #
  verts = getVariable.(dfg, xsyms)
  return doautoinit!(dfg, verts, solveKey=solveKey, singles=singles, N=N, logger=logger)
end
function doautoinit!( dfg::T,
                      xsym::Symbol;
                      solveKey::Symbol=:default,
                      singles::Bool=true,
                      N::Int=100,
                      logger=ConsoleLogger()  )::Bool where T <: AbstractDFG
  #
  return doautoinit!(dfg, [getVariable(dfg, xsym);], solveKey=solveKey, singles=singles, N=N, logger=logger)
end

"""
    $(TYPEDSIGNATURES)

Method to manually initialize a variable using a set of points.

Notes
- Disable automated graphinit on `addFactor!(fg, ...; graphinit=false)
  - any un-initialized variables will automatically be initialized by `solveTree!`

Example:

```julia
# some variable is added to fg
addVariable!(fg, :somepoint3, ContinuousEuclid{2})

# data is organized as (row,col) == (dimension, samples)
pts = randn(2,100)
initManual!(fg, :somepoint3, pts)

# manifold management should be done automatically.
# note upgrades are coming to consolidate with Manifolds.jl, see RoME #244

## it is also possible to initManual! by using existing factors, e.g.
initManual!(fg, :x3, [:x2x3f1])
```

DevNotes
- TODO better document graphinit and treeinit.
"""
function initManual!( variable::DFGVariable, 
                      ptsArr::ManifoldKernelDensity,
                      solveKey::Symbol=:default;
                      dontmargin::Bool=false,
                      N::Int=100 )
  #
  @debug "initManual! $label"
  if !(solveKey in listSolveKeys(variable))
    @debug "$(getLabel(variable)) needs new VND solveKey=$(solveKey)"
    varType = getVariableType(variable)
    setDefaultNodeData!(variable, 0, N, getDimension(varType), solveKey=solveKey, 
    initialized=false, varType=varType, dontmargin=dontmargin)
  end
  setValKDE!(variable, ptsArr, true, solveKey=solveKey)
  return nothing
end
function initManual!( dfg::AbstractDFG, 
                      label::Symbol, 
                      belief::ManifoldKernelDensity,
                      solveKey::Symbol=:default;
                      dontmargin::Bool=false,
                      N::Int=getSolverParams(dfg).N  )
  #
  variable = getVariable(dfg, label)
  initManual!(variable, belief, solveKey, dontmargin=dontmargin, N=N)
  return nothing
end
function initManual!( dfg::AbstractDFG, 
                      label::Symbol, 
                      usefcts::Vector{Symbol},
                      solveKey::Symbol=:default;
                      dontmargin::Bool=false,
                      N::Int=getSolverParams(dfg).N )
  #
  pts = predictbelief(dfg, label, usefcts, solveKey=solveKey)[1]
  vert = getVariable(dfg, label)
  Xpre = AMP.manikde!(pts, getVariableType(vert) |> getManifolds )
  initManual!(vert, Xpre, solveKey, dontmargin=dontmargin, N=N )
  # setValKDE!(vert, Xpre, true, solveKey=solveKey)
  # return nothing
end


function initManual!( dfg::AbstractDFG, 
                      sym::Symbol, 
                      pts::AbstractVector{P}, 
                      solveKey::Symbol=:default ) where {P}
  #
  var = getVariable(dfg, sym)
  pp = manikde!(pts, getManifold(var))
  initManual!(var,pp, solveKey)
end

const initVariableManual! = initManual!

"""
    $SIGNATURES

Set solveKey values of `dest::AbstractDFG` according to `initKey::Symbol=:graphinit` values.

Notes
- Some flexibility for using two DFGs and different key values, see Examples and code for details.
- Can also be specific with `varList::Vector{Symbol}`.
- Returns `dest` graph.
- Uses the supersolve mechanism.

Examples
```julia
resetInitialValues!(fg)
resetInitialValues!(fg1,fg2)  # into 1 from 2
resetInitialValues!(fg1,fg1,:myotherinit)  # use different init value into solveKey :default
resetInitialValues!(fg1,fg1,:graphinit, :mysolver) # not into solveKey=:default but :mysolver
resetInitialValues!(fg1, varList=[:x1;:l3])  # Specific variables only

# Into `fgNew` object, leaving `fg` untouched
fgNew = deepcopy(fg)
resetInitialValues!(fgNew,fg)
```

Related

initManual!, graphinit (keyword)
"""
function resetInitialValues!(dest::AbstractDFG,
                             src::AbstractDFG=dest,
                             initKey::Symbol=:graphinit,
                             solveKey::Symbol=:default;
                             varList::AbstractVector{Symbol}=ls(dest))
  #
  for vs in varList
    vnd = getSolverData(getVariable(src, vs), initKey)
    # guess we definitely want to use copy to preserve the initKey memory
    updateVariableSolverData!(dest,vs,vnd,solveKey,true; warn_if_absent=false)
  end
  return dest
end
const resetInitValues! = resetInitialValues!

"""
    $SIGNATURES

Ensure that no variables set as `solvable=1` are floating free without any connected `solvable=1` factors.  If any found, then set those 'free' variable's `solvable=solvableFallback` (default `0`).

Related

ensureAllInitialized!
"""
function ensureSolvable!(dfg::AbstractDFG; solvableTarget::Int=1, solvableFallback::Int=0)
  # workaround in case isolated variables occur
  solvVars = ls(dfg, solvable=solvableTarget)
  varHasFact = (x->length(ls(dfg,x, solvable=solvableTarget))==0).(solvVars)
  blankVars = solvVars[findall(varHasFact)]
  if 0 < length(blankVars)
    @warn("solveTree! dissallows solvable variables without any connected solvable factors -- forcing solvable=0 on $(blankVars)")
    (x->setSolvable!(dfg, x, solvableFallback)).(blankVars)
  end
  return blankVars
end

"""
    $SIGNATURES

Perform `graphinit` over all variables with `solvable=1` (default).

Related

ensureSolvable!, (EXPERIMENTAL 'treeinit')
"""
function initAll!(dfg::AbstractDFG,
                  solveKey::Symbol=:default; 
                  solvable::Int=1,
                  N::Int=getSolverParams(dfg).N )
  #
  # allvarnodes = getVariables(dfg)
  syms = intersect(getAddHistory(dfg), ls(dfg, solvable=solvable) )
  # syms = ls(dfg, solvable=solvable) # |> sortDFG
  
  # May have to first add the solveKey VNDs if they are not yet available
  for sym in syms
    var = getVariable(dfg, sym)
    # does SolverData exist for this solveKey?
    if !( solveKey in listSolveKeys(var) )
      varType = getVariableType(var)
      # accept complete defaults for a novel solveKey
      setDefaultNodeData!(var, 0, N, getDimension(varType), solveKey=solveKey, 
                          initialized=false, varType=varType, dontmargin=false)
    end
  end

  # do the init
  repeatCount = 0
  repeatFlag = true
  while repeatFlag
    repeatFlag = false
    repeatCount += 1
    if 10 < repeatCount
      @info "not able to initialize all variables via the factor graph, abort autoinit."
      break;
    end
    for sym in syms
      var = getVariable(dfg, sym)
      # is this SolverData initialized?
      if !isInitialized(var, solveKey)
        @info "$(var.label) is not initialized, and will do so now..."
        doautoinit!(dfg, [var;], solveKey=solveKey, singles=true)
        !isInitialized(var, solveKey) ? (repeatFlag = true) : nothing
      end
    end
  end
  nothing
end

function assembleFactorName(dfg::AbstractDFG,
                            Xi::Vector{<:DFGVariable} )
  #

  existingFactorLabels = listFactors(dfg)
  existingFactorLabelDict = Dict(existingFactorLabels .=> existingFactorLabels)
  namestring = ""
  for vert in Xi #f.Xi
    namestring = string(namestring,vert.label)
  end
  opt = getSolverParams(dfg)
  for i in 1:opt.maxincidence
    tempnm = string(namestring, "f$i")
    if !haskey(existingFactorLabelDict, Symbol(tempnm))
      namestring = tempnm
      break
    end
    i != opt.maxincidence ? nothing : error("Artificial restriction to not connect more than $(opt.maxincidence) factors to a variable (bad for sparsity), try setting getSolverParams(fg).maxincidence=1000 to adjust this restriction.")
  end
  return Symbol(namestring)
end

"""
    $(SIGNATURES)

Add factor with user defined type <: FunctorInferenceType to the factor graph
object. Define whether the automatic initialization of variables should be
performed.  Use order sensitive `multihypo` keyword argument to define if any
variables are related to data association uncertainty.

Experimental
- `inflation`, to better disperse kernels before convolution solve, see IIF #1051.
"""
function addFactor!(dfg::AbstractDFG,
                    Xi::Vector{<:DFGVariable},
                    usrfnc::R;
                    multihypo::Vector{Float64}=Float64[],
                    nullhypo::Float64=0.0,
                    solvable::Int=1,
                    tags::Vector{Symbol}=Symbol[],
                    timestamp::Union{DateTime,ZonedDateTime}=now(localzone()),
                    graphinit::Bool=getSolverParams(dfg).graphinit,
                    threadmodel=SingleThreaded,
                    suppressChecks::Bool=false,
                    inflation::Real=getSolverParams(dfg).inflation  ) where
                      {R <: FunctorInferenceType}
  #
  # depcrecation

  varOrderLabels = Symbol[v.label for v=Xi]
  namestring = assembleFactorName(dfg, Xi)
  solverData = getDefaultFactorData(dfg, 
                                    Xi, 
                                    deepcopy(usrfnc), 
                                    multihypo=multihypo, 
                                    nullhypo=nullhypo, 
                                    threadmodel=threadmodel,
                                    inflation=inflation)
  newFactor = DFGFactor(Symbol(namestring),
                        varOrderLabels,
                        solverData;
                        tags=Set(union(tags, [:FACTOR])),
                        solvable=solvable,
                        timestamp=timestamp)
  #

  success = DFG.addFactor!(dfg, newFactor)

  # TODO: change this operation to update a conditioning variable
  graphinit && doautoinit!(dfg, Xi, singles=false)

  return newFactor
end

function DFG.addFactor!(dfg::AbstractDFG,
                        xisyms::Vector{Symbol},
                        usrfnc::FunctorInferenceType;
                        multihypo::Vector{<:Real}=Float64[],
                        nullhypo::Float64=0.0,
                        solvable::Int=1,
                        timestamp::Union{DateTime,ZonedDateTime}=now(localzone()),
                        tags::Vector{Symbol}=Symbol[],
                        graphinit::Bool=getSolverParams(dfg).graphinit,
                        threadmodel=SingleThreaded,
                        inflation::Real=getSolverParams(dfg).inflation,
                        suppressChecks::Bool=false  )
  #
  # depcrecation

  # basic sanity check for unary vs n-ary
  if !suppressChecks && length(xisyms) == 1 && !(usrfnc isa AbstractPrior) && !(usrfnc isa Mixture)
    @warn("Listing only one variable $xisyms for non-unary factor type $(typeof(usrfnc))")
  end

  variables = getVariable.(dfg, xisyms)
  # verts = map(vid -> DFG.getVariable(dfg, vid), xisyms)
  addFactor!(dfg, variables, usrfnc, multihypo=multihypo, nullhypo=nullhypo, solvable=solvable, tags=tags, graphinit=graphinit, threadmodel=threadmodel, timestamp=timestamp, inflation=inflation )
end




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
- default is to use `solvable=1` and ignore factors and variables that might be used for dead reckoning or similar.

Future
- TODO: `A` should be sparse data structure (when we exceed 10'000 var dims)
- TODO: Incidence matrix is rectagular and adjacency is the square.
"""
function getEliminationOrder( dfg::G;
                              ordering::Symbol=:qr,
                              solvable::Int=1,
                              constraints::Vector{Symbol}=Symbol[]) where G <: AbstractDFG
  #
  @assert 0 == length(constraints) || ordering == :ccolamd "Must use ordering=:ccolamd when trying to use constraints"
  # Get the sparse adjacency matrix, variable, and factor labels
  adjMat, permuteds, permutedsf = DFG.getBiadjacencyMatrix(dfg, solvable=solvable)
  # adjMat, permuteds, permutedsf = DFG.getAdjacencyMatrixSparse(dfg, solvable=solvable)

  # Create dense adjacency matrix

  p = Int[]
  if ordering==:chol
    # hack for dense matrix....
    A = adjMat
    p = cholesky(Matrix(A'A),Val(true)).piv
    @warn "check that cholesky ordering is not reversed -- basically how much fill in (separator size) are you seeing???  Long skinny chains in tree is bad."
  elseif ordering==:qr
    # hack for dense matrix....
    A = Array(adjMat)
    # this is the default
    q,r,p = qr(A, Val(true))
    p .= p |> reverse
  elseif ordering==:ccolamd
    cons = zeros(SuiteSparse_long, length(adjMat.colptr) - 1)
    cons[findall(x->x in constraints, permuteds)] .= 1
    p = Ccolamd.ccolamd(adjMat, cons)
    @warn "Ccolamd is experimental in IIF at this point in time."
  else
    prtslperr("getEliminationOrder -- cannot do the requested ordering $(ordering)")
  end

  # Return the variable ordering that we should use for the Bayes map
  # reverse order checked in #475 and #499
  return permuteds[p]
end


# lets create all the vertices first and then deal with the elimination variables thereafter
function addBayesNetVerts!( dfg::AbstractDFG,
                            elimOrder::Array{Symbol,1} )
  #
  for pId in elimOrder
    vert = DFG.getVariable(dfg, pId)
    if getSolverData(vert).BayesNetVertID == nothing || getSolverData(vert).BayesNetVertID == :_null # Special serialization case of nothing
      @debug "[AddBayesNetVerts] Assigning $pId.data.BayesNetVertID = $pId"
      getSolverData(vert).BayesNetVertID = pId
    else
      @warn "addBayesNetVerts -- Something is wrong, variable '$pId' should not have an existing Bayes net reference to '$(getSolverData(vert).BayesNetVertID)'"
    end
  end
end

function addConditional!( dfg::AbstractDFG,
                          vertId::Symbol,
                          Si::Vector{Symbol} )
  #
  bnv = DFG.getVariable(dfg, vertId)
  bnvd = getSolverData(bnv)
  bnvd.separator = Si
  for s in Si
    push!(bnvd.BayesNetOutVertIDs, s)
  end
  return nothing
end

function addChainRuleMarginal!( dfg::AbstractDFG,
                                Si::Vector{Symbol} )
  #

  lbls = String[]
  genmarg = GenericMarginal()
  Xi = map(v -> DFG.getVariable(dfg, v), Si)
  # @info "adding marginal to"
  # for x in Xi
  #   @info "x.index=",x.index
  # end
  addFactor!( dfg, Xi, genmarg, graphinit=false, suppressChecks=true )
  nothing
end

function rmVarFromMarg( dfg::AbstractDFG,
                        fromvert::DFGVariable,
                        gm::Vector{DFGFactor}  )
  #

  @debug " - Removing $(fromvert.label)"
  for m in gm
    @debug "Looking at $(m.label)"
    for n in DFG.getNeighbors(dfg, m) #x1, x2
      if n == fromvert.label # n.label ==? x1
        @debug "   - Breaking link $(m.label)->$(fromvert.label)..."
        @debug "     - Original links: $(DFG.ls(dfg, m))"
        remvars = setdiff(DFG.ls(dfg, m), [fromvert.label])
        @debug "     - New links: $remvars"

        DFG.deleteFactor!(dfg, m) # Remove it
        if length(remvars) > 0
          @debug "$(m.label) still has links to other variables, readding it back..."
          addFactor!(dfg, remvars, _getCCW(m).usrfnc!, graphinit=false, suppressChecks=true )
        else
          @debug "$(m.label) doesn't have any other links, not adding it back..."
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

function buildBayesNet!(dfg::AbstractDFG,
                        elimorder::Vector{Symbol};
                        solvable::Int=1 )
  #
  # addBayesNetVerts!(dfg, elimorder)
  for v in elimorder
    @debug """ 
                Eliminating $(v)
                ===============
          """
    # which variable are we eliminating

    # all factors adjacent to this variable
    fi = Symbol[]
    Si = Symbol[]
    gm = DFGFactor[]

    vert = DFG.getVariable(dfg, v)
    for fctId in DFG.getNeighbors(dfg, vert, solvable=solvable)
      fct = DFG.getFactor(dfg, fctId)
      if (getSolverData(fct).eliminated != true)
        push!(fi, fctId)
        for sepNode in DFG.getNeighbors(dfg, fct, solvable=solvable)
          # TODO -- validate !(sepNode.index in Si) vs. older !(sepNode in Si)
          if sepNode != v && !(sepNode in Si) # Symbol comparison!
            push!(Si,sepNode)
          end
        end
        getSolverData(fct).eliminated = true
      end

      if typeof(_getCCW(fct)) == CommonConvWrapper{GenericMarginal}
        push!(gm, fct)
      end
    end

    if v != elimorder[end]
      addConditional!(dfg, v, Si)
      # not yet inserting the new prior p(Si) back into the factor graph
    end

    # mark variable
    getSolverData(vert).eliminated = true

    # TODO -- remove links from current vertex to any marginals
    rmVarFromMarg(dfg, vert, gm)

    #add marginal on remaining variables... ? f(xyz) = f(x | yz) f(yz)
    # new function between all Si (round the outside, right the outside)
    length(Si) > 0 && addChainRuleMarginal!(dfg, Si)

  end
  return nothing
end



"""
    $(SIGNATURES)

Get KernelDensityEstimate kde estimate stored in variable node.
"""
getBelief(vnd::VariableNodeData) = manikde!(getVal(vnd), getBW(vnd)[:,1], getVariableType(vnd) )

getBelief(v::DFGVariable, solvekey::Symbol=:default) = getBelief(getSolverData(v, solvekey))
getBelief(dfg::AbstractDFG, lbl::Symbol, solvekey::Symbol=:default) = getBelief(getVariable(dfg, lbl), solvekey)

const getKDE = getBelief


#
