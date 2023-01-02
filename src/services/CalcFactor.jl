# New factor interface, something perhaps like this

export calcFactorResidualTemporary

# NOTE, the full concrete type is recovered in reconstFactorData
getFactorOperationalMemoryType(dfg::SolverParams) = CommonConvWrapper
# difficult type piracy case needing both types NoSolverParams and CommonConvWrapper.
getFactorOperationalMemoryType(dfg::NoSolverParams) = CommonConvWrapper

getManifold(fct::DFGFactor{<:CommonConvWrapper}) = getManifold(_getCCW(fct))

function _getDimensionsPartial(ccw::CommonConvWrapper)
  # @warn "_getDimensionsPartial not ready for use yet"
  return ccw.partialDims
end
function _getDimensionsPartial(data::GenericFunctionNodeData)
  return _getCCW(data) |> _getDimensionsPartial
end
_getDimensionsPartial(fct::DFGFactor) = _getDimensionsPartial(_getCCW(fct))
function _getDimensionsPartial(fg::AbstractDFG, lbl::Symbol)
  return _getDimensionsPartial(getFactor(fg, lbl))
end

# Helper function to construct CF from a CCW
function CalcFactor(
  ccwl::CommonConvWrapper;
  factor = ccwl.usrfnc!,
  _sampleIdx = 0,
  _legacyParams = ccwl.varValsAll,
  _allowThreads = true,
  cache = ccwl.dummyCache,
  fullvariables = ccwl.fullvariables,
  solvefor = ccwl.varidx,
  manifold = getManifold(ccwl)
)
  #
  # FIXME using ccwl.dummyCache is not thread-safe
  return CalcFactor(
    factor,
    _sampleIdx,
    _legacyParams,
    _allowThreads,
    cache,
    tuple(fullvariables...),
    solvefor,
    manifold
  )
end

"""
    $SIGNATURES

Sample the factor stochastic model `N::Int` times and store the samples in the preallocated `ccw.measurement` container.

DevNotes
- Use in place operations where possible and remember `measurement` is a `::Tuple`.
- TODO only works on `.threadid()==1` at present, see #1094
- Also see, JuliaRobotics/RoME.jl#465
"""
sampleFactor(cf::CalcFactor{<:AbstractFactor}, N::Int = 1) = [getSample(cf) for _ = 1:N]

function Base.show(io::IO, x::CalcFactor)
  println(io)
  printstyled(io, " CalcFactor:\n"; color = :blue)
  return println(io, "  .factor: ", typeof(x.factor))
end

Base.show(io::IO, ::MIME"text/plain", x::CalcFactor) = show(io, x)

"""
    $SIGNATURES

Function to calculate measurement dimension from factor sampling.

Notes
- Will not work in all situations, but good enough so far.
  - # TODO standardize via domain or manifold definition...??
"""
function calcZDim(cf::CalcFactor{T}) where {T <: AbstractFactor}
  #
  M = getManifold(cf) # getManifold(T)
  try
    return manifold_dimension(M)
  catch
    @warn "no method getManifold(::$(string(T))), calcZDim will attempt legacy length(sample) method instead"
  end

  # NOTE try to make sure we get matrix back (not a vector)
  smpls = sampleFactor(cf, 2)[1]
  return length(smpls)
end

calcZDim(ccw::CommonConvWrapper) = calcZDim(CalcFactor(ccw))

calcZDim(cf::CalcFactor{<:GenericMarginal}) = 0

calcZDim(cf::CalcFactor{<:ManifoldPrior}) = manifold_dimension(cf.manifold)

"""
    $SIGNATURES

Helper function for evaluating factor residual functions, by adding necessary `CalcFactor` wrapper.
  
Notes
- Factor must already be in a factor graph to work
- Will not yet properly support all multihypo nuances, more a function for testing
- Useful for debugging a factor. 

Example
```julia
fg = generateGraph_Kaess()

residual = calcFactorResidual(fg, :x1x2f1, [1.0], [0.0], [0.0])
```

Related

[`calcFactorResidualTemporary`](@ref), [`_evalFactorTemporary!`](@ref), [`evalFactor`](@ref), [`approxConv`](@ref)
"""
function calcFactorResidual(
  dfgfct::DFGFactor,
  args...;
  ccw::CommonConvWrapper = IIF._getCCW(dfgfct),
)
  return CalcFactor(ccw)(args...)
end
function calcFactorResidual(dfg::AbstractDFG, fctsym::Symbol, args...)
  return calcFactorResidual(getFactor(dfg, fctsym), args...)
end

"""
    $SIGNATURES

Evaluate the residual function for a single sample.

Notes
- Binary factors only at this stage, and `multihypo` does not have to be considered in this test
- Assumes calculation is for a single particle, so `meas::Tuple{Z,other}` is only a single particles value.

Example
```julia
residual = calcFactorResidualTemporary(Pose2Pose2(...), (RoME.Pose2,RoME.Pose2), (z_i,), (x1, x2))
```

Related

[`calcFactorResidual`](@ref), [`CalcResidual`](@ref), [`_evalFactorTemporary!`](@ref), [`approxConv`](@ref), [`_buildGraphByFactorAndTypes!`](@ref)
"""
function calcFactorResidualTemporary(
  fct::AbstractRelative,
  varTypes::Tuple,
  measurement,
  pts::Tuple;
  tfg::AbstractDFG = initfg(),
  _blockRecursion::Bool = false,
  doTime::Bool = false,
)
  #

  # build a new temporary graph
  _, _dfgfct = _buildGraphByFactorAndTypes!(
    fct,
    varTypes,
    pts;
    dfg = tfg,
    _blockRecursion = _blockRecursion,
  )

  # get a fresh measurement if needed
  _measurement = if measurement != [] #length(measurement) != 0
    measurement
  else
    # now use the CommonConvWrapper object in `_dfgfct`
    cfo = CalcFactor(_getCCW(_dfgfct))
    sampleFactor(cfo, 1)[1]
  end

  # assume a single sample point is being run
  res = if doTime
    @time res = calcFactorResidual(_dfgfct, _measurement, pts...)
    res
  else
    calcFactorResidual(_dfgfct, _measurement, pts...)
  end
  return res
end

## =============================================================================================
## FactorOperationalMemory helper constructors
## =============================================================================================

function CommonConvWrapper(
  usrfnc::T,
  fullvariables, #::Tuple ::Vector{<:DFGVariable};
  varValsAll::Tuple,
  X::AbstractVector{P}; #TODO remove X completely
  # xDim::Int = size(X, 1),
  userCache::CT = nothing,
  manifold = getManifold(usrfnc),
  partialDims::AbstractVector{<:Integer} = 1:length(X),
  partial::Bool = false,
  nullhypo::Real = 0,
  inflation::Real = 3.0,
  hypotheses::H = nothing,
  certainhypo = nothing,
  activehypo = collect(1:length(varValsAll)),
  measurement::AbstractVector = Vector(Vector{Float64}()),
  varidx::Int = 1,
  particleidx::Int = 1,
  res::AbstractVector{<:Real} = zeros(manifold_dimension(manifold)), # zDim
  gradients = nothing,
) where {T <: AbstractFactor, P, H, CT}
  #
  return CommonConvWrapper(
    usrfnc,
    tuple(fullvariables...),
    varValsAll,
    userCache,
    manifold,
    partialDims,
    partial,
    # xDim,
    Float64(nullhypo),
    inflation,
    hypotheses,
    certainhypo,
    activehypo,
    measurement,
    varidx,
    particleidx,
    res,
    gradients,
  )
end

# the same as legacy, getManifold(ccwl.usrfnc!)
getManifold(ccwl::CommonConvWrapper) = ccwl.manifold
getManifold(cf::CalcFactor) = cf.manifold

function _resizePointsVector!(
  vecP::AbstractVector{P},
  mkd::ManifoldKernelDensity,
  N::Int,
) where {P}
  #
  pN = length(vecP)
  resize!(vecP, N)
  for j = pN:N
    smp = AMP.sample(mkd, 1)[1]
    # @show j, smp, typeof(smp), typeof(vecP[j])
    vecP[j] = smp[1]
  end

  return vecP
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

DevNotes
- FIXME ARR internally should become a NamedTuple
"""
function _prepParamVec(
  Xi::Vector{<:DFGVariable},
  solvefor::Union{Nothing, Symbol},
  N::Int = 0;
  solveKey::Symbol = :default,
)
  #
  # FIXME refactor to new NamedTuple instead
  varParamsAll = getVal.(Xi; solveKey)
  # Xi_labels = getLabel.(Xi)
  sfidx = if isnothing(solvefor)
    0
  else
    findfirst(==(solvefor), getLabel.(Xi))
  end
  # sfidx = isnothing(sfidx) ? 0 : sfidx

  # this line does nothing...
  # maxlen = N # FIXME see #105

  LEN = length.(varParamsAll)
  maxlen = maximum([N; LEN])

  # resample variables with too few kernels (manifolds points)
  SAMP = LEN .< maxlen
  for i = 1:length(Xi)
    if SAMP[i]
      Pr = getBelief(Xi[i], solveKey)
      _resizePointsVector!(varParamsAll[i], Pr, maxlen)
    end
  end

  # TODO --rather define reusable memory for the proposal
  # we are generating a proposal distribution, not direct replacement for existing memory and hence the deepcopy.
  if sfidx > 0
    varParamsAll[sfidx] = deepcopy(varParamsAll[sfidx])
  end

  varValsAll = tuple(varParamsAll...)
  # FIXME, forcing maxlen to N results in errors (see test/testVariousNSolveSize.jl) see #105
  # maxlen = N == 0 ? maxlen : N
  return varValsAll, maxlen, sfidx
end

"""
    $SIGNATURES
Internal method to set which dimensions should be used as the decision variables for later numerical optimization.
"""
function _setCCWDecisionDimsConv!(
  ccwl::Union{CommonConvWrapper{F}, CommonConvWrapper{Mixture{N_, F, S, T}}},
  xDim::Int
) where {
  N_,
  F <: Union{
    AbstractManifoldMinimize,
    AbstractRelativeMinimize,
    AbstractRelativeRoots,
    AbstractPrior,
  },
  S,
  T,
}
  #

  # NOTE should only be done in the constructor
  ccwl.partialDims = if ccwl.partial
    Int[ccwl.usrfnc!.partial...]
  else
    # NOTE this is the target variable dimension (not factor manifold dimension) 
    Int[1:xDim...] # ccwl.xDim
  end

  return nothing
end

function attemptGradientPrep(
  varTypes,
  usrfnc,
  varParamsAll,
  multihypo,
  meas_single,
  _blockRecursion,
)
  # prepare new cached gradient lambdas (attempt)
  try
    # https://github.com/JuliaRobotics/IncrementalInference.jl/blob/db7ff84225cc848c325e57b5fb9d0d85cb6c79b8/src/DispatchPackedConversions.jl#L46
    # also https://github.com/JuliaRobotics/DistributedFactorGraphs.jl/issues/590#issuecomment-891450762
    # FIXME, suppressing nested gradient propagation on GenericMarginals for the time being, see #1010
    if (!_blockRecursion) && usrfnc isa AbstractRelative && !(usrfnc isa GenericMarginal)
      # take first value from each measurement-tuple-element
      measurement_ = meas_single
      # compensate if no info available during deserialization
      # take the first value from each variable param
      pts_ = map(x -> x[1], varParamsAll)
      # FIXME, only using first meas and params values at this time...
      # NOTE, must block recurions here, since FGC uses this function to calculate numerical gradients on a temp fg.
      # assume for now fractional-var in multihypo have same varType
      hypoidxs = _selectHypoVariables(pts_, multihypo)
      gradients = FactorGradientsCached!(
        usrfnc,
        tuple(varTypes[hypoidxs]...),
        measurement_,
        tuple(pts_[hypoidxs]...);
        _blockRecursion = true,
      )

      return gradients
    end
  catch e
    @warn "Unable to create measurements and gradients for $usrfnc during prep of CCW, falling back on no-partial information assumption.  Enable ENV[\"JULIA_DEBUG\"] = \"IncrementalInference\" for @debug printing to see the error."
    # rethrow(e)
    @debug(e)
  end
  return nothing
end

"""
    $SIGNATURES

Notes
- Can be called with `length(Xi)==0`
"""
function _prepCCW(
  Xi::Vector{<:DFGVariable},
  usrfnc::T;
  multihypo::Union{Nothing, <:Distributions.Categorical} = nothing,
  nullhypo::Real = 0.0,
  certainhypo = if multihypo !== nothing
    collect(1:length(multihypo.p))[multihypo.p .== 0.0]
  else
    collect(1:length(Xi))
  end,
  inflation::Real = 0.0,
  solveKey::Symbol = :default,
  _blockRecursion::Bool = false,
  userCache::CT = nothing,
) where {T <: AbstractFactor, CT}
  #
  if length(Xi) !== 0
    nothing
  else
    @debug("cannot prep ccw.param list with length(Xi)==0, see DFG #590")
  end

  # TODO check no Anys, see #1321
  _varValsAll, maxlen, sfidx = _prepParamVec(Xi, nothing, 0; solveKey) # Nothing for init.

  manifold = getManifold(usrfnc)
  # standard factor metadata
  solvefor = length(Xi)
  fullvariables = tuple(Xi...) # convert(Vector{DFGVariable}, Xi)
  # create a temporary CalcFactor object for extracting the first sample
  # TODO, deprecate this:  guess measurement points type
  # MeasType = Vector{Float64} # FIXME use `usrfnc` to get this information instead
  _cf = CalcFactor(
    usrfnc,
    0,
    _varValsAll,
    false,
    userCache,
    fullvariables,
    solvefor,
    manifold
  )

  # get a measurement sample
  meas_single = sampleFactor(_cf, 1)[1]

  elT = typeof(meas_single)

  #TODO preallocate measurement?
  measurement = Vector{elT}()

  # partialDims are sensitive to both which solvefor variable index and whether the factor is partial
  partial = hasfield(T, :partial) # FIXME, use isPartial function instead
  partialDims = if partial
    Int[usrfnc.partial...]
  else
    Int[]
  end

  # FIXME, should incorporate multihypo selection
  varTypes = getVariableType.(fullvariables)

  # as per struct CommonConvWrapper
  gradients = attemptGradientPrep(
    varTypes,
    usrfnc,
    _varValsAll,
    multihypo,
    meas_single,
    _blockRecursion,
  )

  # variable Types
  pttypes = getVariableType.(Xi) .|> getPointType
  PointType = 0 < length(pttypes) ? pttypes[1] : Vector{Float64}
  if !isconcretetype(PointType)
    @warn "_prepCCW PointType is not concrete $PointType" maxlog=50
  end

  return CommonConvWrapper(
    usrfnc,
    fullvariables,
    _varValsAll,
    PointType[];
    userCache, # should be higher in args list
    manifold,  # should be higher in args list
    partialDims,
    partial,
    nullhypo,
    inflation,
    hypotheses = multihypo,
    certainhypo,
    measurement,
    gradients,
  )
end

function updateMeasurement!(
  ccwl::CommonConvWrapper,
  N::Int=1;
  measurement::AbstractVector = Vector{Tuple{}}(),
  needFreshMeasurements::Bool=true,
  _allowThreads::Bool = true
)
  # FIXME do not divert Mixture for sampling
  
  # option to disable fresh samples or user provided
  if needFreshMeasurements
    # TODO this is only one thread, make this a for loop for multithreaded sampling
    sampleFactor!(ccwl, N; _allowThreads)
  elseif 0 < length(measurement) 
    resize!(ccwl.measurement, length(measurement))
    ccwl.measurement[:] = measurement
  end

  nothing
end

"""
    $(SIGNATURES)

Prepare a common functor computation object `prepareCommonConvWrapper{T}` containing 
the user factor functor along with additional variables and information using during 
approximate convolution computations.

DevNotes
- TODO consolidate with others, see https://github.com/JuliaRobotics/IncrementalInference.jl/projects/6
"""
function _updateCCW!(
  F_::Type{<:AbstractRelative},
  ccwl::CommonConvWrapper{F},
  Xi::AbstractVector{<:DFGVariable},
  solvefor::Symbol,
  N::Integer;
  measurement = Vector{Tuple{}}(),
  needFreshMeasurements::Bool = true,
  solveKey::Symbol = :default,
) where {F <: AbstractFactor} # F might be Mixture
  #
  if length(Xi) !== 0
    nothing
  else
    @debug("cannot prep ccw.param list with length(Xi)==0, see DFG #590")
  end

  # FIXME, order of fmd ccwl cf are a little weird and should be revised.
  # FIXME maxlen should parrot N (barring multi-/nullhypo issues)
  _varValsAll, maxlen, sfidx = _prepParamVec(Xi, solvefor, N; solveKey)

  # TODO, ensure all values (not just multihypothesis) is correctly used from here
  for (i,varVal) in enumerate(_varValsAll)
    resize!(ccwl.varValsAll[i],length(varVal))
    ccwl.varValsAll[i][:] = varVal
  end

  # set the 'solvefor' variable index -- i.e. which connected variable of the factor is being computed in this convolution. 
  ccwl.varidx = sfidx

  # TODO better consolidation still possible
  # FIXME ON FIRE, what happens if this is a partial dimension factor?  See #1246
  # FIXME, confirm this is hypo sensitive selection from Xi, better to use double indexing for clarity getDimension(ccw.fullvariables[hypoidx[sfidx]])
  xDim = getDimension(getVariableType(Xi[sfidx])) # ccwl.varidx
  # ccwl.xDim = xDim
  # TODO maybe refactor different type or api call?

  # setup the partial or complete decision variable dimensions for this ccwl object
  # NOTE perhaps deconv has changed the decision variable list, so placed here during consolidation phase
  # TODO, should this not be part of `prepareCommonConvWrapper` -- only here do we look for .partial
  _setCCWDecisionDimsConv!(ccwl, xDim)

  # FIXME do not divert Mixture for sampling

  updateMeasurement!(ccwl, maxlen; needFreshMeasurements, measurement, _allowThreads=true)

  # used in ccw functor for AbstractRelativeMinimize
  resize!(ccwl.res, _getZDim(ccwl))
  fill!(ccwl.res, 0.0)

  # calculate new gradients
  # J = ccwl.gradients(measurement..., pts...)

  return sfidx, maxlen
end

function _updateCCW!(
  F_::Type{<:AbstractPrior},
  ccwl::CommonConvWrapper{F},
  Xi::AbstractVector{<:DFGVariable},
  solvefor::Symbol,
  N::Integer;
  measurement = Vector{Tuple{}}(),
  needFreshMeasurements::Bool = true,
  solveKey::Symbol = :default,
) where {F <: AbstractFactor} # F might be Mixture
  # FIXME, NEEDS TO BE CLEANED UP AND WORK ON MANIFOLDS PROPER
  # fnc = ccwl.usrfnc!
  sfidx = findfirst(getLabel.(Xi) .== solvefor)
  @assert sfidx == 1 "Solving on Prior with CCW should have sfidx=1, priors are unary factors."
  # sfidx = 1 #  why hardcoded to 1, maybe for only the AbstractPrior case here

  # setup the partial or complete decision variable dimensions for this ccwl object
  # NOTE perhaps deconv has changed the decision variable list, so placed here during consolidation phase
  _setCCWDecisionDimsConv!(ccwl, getDimension(getVariableType(Xi[sfidx])))

  solveForPts = getVal(Xi[sfidx]; solveKey)
  maxlen = maximum([N; length(solveForPts); length(ccwl.varValsAll[sfidx])])  # calcZDim(ccwl); length(measurement[1])

  # FIXME do not divert Mixture for sampling
  # update ccwl.measurement values
  updateMeasurement!(ccwl, maxlen; needFreshMeasurements, measurement, _allowThreads=true)

  return sfidx, maxlen
end

# TODO, can likely deprecate this
function _updateCCW!(
  ccwl::Union{CommonConvWrapper{F}, CommonConvWrapper{Mixture{N_, F, S, T}}},
  Xi::AbstractVector{<:DFGVariable},
  solvefor::Symbol,
  N::Integer;
  kw...,
) where {N_, F <: AbstractRelative, S, T}
  #
  return _updateCCW!(F, ccwl, Xi, solvefor, N; kw...)
end

function _updateCCW!(
  ccwl::Union{CommonConvWrapper{F}, CommonConvWrapper{Mixture{N_, F, S, T}}},
  Xi::AbstractVector{<:DFGVariable},
  solvefor::Symbol,
  N::Integer;
  kw...,
) where {N_, F <: AbstractPrior, S, T}
  #
  return _updateCCW!(F, ccwl, Xi, solvefor, N; kw...)
end


#
