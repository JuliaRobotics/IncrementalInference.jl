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
function CalcFactorNormSq(
  ccwl::CommonConvWrapper;
  factor = ccwl.usrfnc!,
  _sampleIdx = ccwl.particleidx[],
  _legacyParams = ccwl.varValsAll[],
  _allowThreads = true,
  cache = ccwl.dummyCache,
  fullvariables = ccwl.fullvariables,
  solvefor = ccwl.varidx[],
  manifold = getManifold(ccwl),
  slack=nothing,
)
  #
  # FIXME using ccwl.dummyCache is not thread-safe
  return CalcFactorNormSq(
    factor,
    _sampleIdx,
    _legacyParams,
    _allowThreads,
    cache,
    tuple(fullvariables...),
    solvefor,
    manifold,
    ccwl.measurement,
    slack
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

calcZDim(ccw::CommonConvWrapper) = calcZDim(CalcFactorNormSq(ccw))

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

[`calcFactorResidualTemporary`](@ref), [`_evalFactorTemporary!`](@ref), [`evalFactor`](@ref), [`approxConvBelief`](@ref)
"""
function calcFactorResidual(
  dfgfct::DFGFactor,
  args...;
  ccw::CommonConvWrapper = IIF._getCCW(dfgfct),
)
  return CalcFactorNormSq(ccw)(args...)
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

[`calcFactorResidual`](@ref), [`CalcResidual`](@ref), [`_evalFactorTemporary!`](@ref), [`approxConvBelief`](@ref), [`_buildGraphByFactorAndTypes!`](@ref)
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
    cfo = CalcFactorNormSq(_getCCW(_dfgfct))
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


# the same as legacy, getManifold(ccwl.usrfnc!)
getManifold(ccwl::CommonConvWrapper) = ccwl.manifold
getManifold(cf::CalcFactor) = getManifold(cf.factor)

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

function _checkVarValPointers(dfg::AbstractDFG, fclb::Symbol)
  vars = getVariable.(dfg, getVariableOrder(dfg,fclb))
  ptrsV = pointer.(getVal.(vars))
  ccw = _getCCW(dfg, fclb)
  ptrsC = pointer.(ccw.varValsAll[])
  ptrsV, ptrsC
end

"""
    $(SIGNATURES)

Prepare the particle arrays `ARR` to be used for approximate convolution.
This function ensures that ARR has te same dimensions among all the parameters.
Function returns with ARR[sfidx] pointing at newly allocated deepcopy of the
existing values in getVal(Xi[.label==solvefor]).

Notes
- 2023Q2 intended use, only create VarValsAll the first time a factor added/reconstructed
- Return values `sfidx` is the element in ARR where `Xi.label==solvefor` and
- `maxlen` is length of all (possibly resampled) `ARR` contained particles.
- `Xi` is order sensitive.
- for initialization, solveFor = Nothing.
- `P = getPointType(<:InferenceVariable)`
"""
function _createVarValsAll(
  variables::AbstractVector{<:DFGVariable};
  solveKey::Symbol = :default,
)
  #
  # Note, NamedTuple once upon a time created way too much recompile load on repeat solves, #1564
  # FIXME ON FIRE issue on deserialization
  valsAll = []

  # when deserializing a factor, a new ccw gets created but the variables may not yet have VND entries
  for var_i in variables
    push!(
      valsAll,
      if haskey(getSolverDataDict(var_i), solveKey)
        getVal(var_i; solveKey)
      else
        Vector{typeof(getPointDefault(getVariableType(var_i)))}()
      end
    )
  end

  varValsAll = tuple(valsAll...)

  # how many points
  LEN = length.(varValsAll)
  maxlen = maximum(LEN)
  # NOTE, forcing maxlen to N results in errors (see test/testVariousNSolveSize.jl) see #105
  # maxlen = N == 0 ? maxlen : N

  # NOTE resize! moves the pointer!!!!!!
    # # allow each variable to have a different number of points, which is resized during compute here
    # # resample variables with too few kernels (manifolds points)
  # SAMP = LEN .< maxlen
  # for i = 1:length(variables)
  #   if SAMP[i]
  #     Pr = getBelief(variables[i], solveKey)
  #     _resizePointsVector!(varValsAll[i], Pr, maxlen)
  #   end
  # end

  # TODO --rather define reusable memory for the proposal
  # we are generating a proposal distribution, not direct replacement for existing memory and hence the deepcopy.
  #  POSSIBLE SOURCE OF HUGE MEMORY CONSUMPTION ALLOCATION

  return varValsAll
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
    AbstractPrior,
  },
  S,
  T,
}
  #

  # NOTE should only be done in the constructor
  newval = if ccwl.partial
    Int[ccwl.usrfnc!.partial...]
  else
    # NOTE this is the target variable dimension (not factor manifold dimension) 
    Int[1:xDim...] # ccwl.xDim
  end
  resize!(ccwl.partialDims, length(newval))
  ccwl.partialDims[:] = newval

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
- _createCCW is likely only used when adding or reconstructing a new factor in the graph,
  - else use _updateCCW 
- Can be called with `length(Xi)==0`
"""
function _createCCW(
  Xi::AbstractVector{<:DFGVariable},
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
  attemptGradients::Bool = true,
  userCache::CT = nothing,
) where {T <: AbstractFactor, CT}
  #
  if length(Xi) !== 0
    nothing
  else
    @debug("cannot prep ccw.param list with length(Xi)==0, see DFG #590")
  end

  # TODO check no Anys, see #1321
  # NOTE, _varValsAll is only a reference to the actual VND.val memory of each variable
  _varValsAll = _createVarValsAll(Xi; solveKey)

  manifold = getManifold(usrfnc)
  # standard factor metadata
  solvefor = length(Xi)
  fullvariables = tuple(Xi...) # convert(Vector{DFGVariable}, Xi)
  # create a temporary CalcFactor object for extracting the first sample

  _cf = CalcFactorNormSq(
    usrfnc,
    1,
    _varValsAll,
    false,
    userCache,
    fullvariables,
    solvefor,
    manifold,
    nothing,
    nothing,
  )

  # get a measurement sample
  meas_single = sampleFactor(_cf, 1)[1]
  elT = typeof(meas_single)
  #TODO preallocate measurement?
  measurement = Vector{elT}()

  #FIXME chicken and egg problem for getting measurement type, so creating twice.
  _cf = CalcFactorNormSq(
    usrfnc,
    1,
    _varValsAll,
    false,
    userCache,
    fullvariables,
    solvefor,
    manifold,
    measurement,
    nothing,
  )


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
  _gradients = if attemptGradients
    attemptGradientPrep(
      varTypes,
      usrfnc,
      _varValsAll,
      multihypo,
      meas_single,
      _blockRecursion,
    )
  else
    nothing
  end

  # variable Types
  pttypes = getVariableType.(Xi) .|> getPointType
  PointType = 0 < length(pttypes) ? pttypes[1] : Vector{Float64}
  if !isconcretetype(PointType)
    @warn "_createCCW PointType is not concrete $PointType" maxlog=50
  end

  # PointType[],
  return CommonConvWrapper(;
    usrfnc! = usrfnc,
    fullvariables,
    varValsAll = Ref(_varValsAll),
    dummyCache = userCache,
    manifold,
    partialDims,
    partial,
    nullhypo = float(nullhypo),
    inflation = float(inflation),
    hyporecipe = HypoRecipeCompute(;
      hypotheses = multihypo,
      certainhypo,
    ),
    measurement,
    _gradients,
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
function _beforeSolveCCW!(
  F_::Type{<:AbstractRelative},
  ccwl::CommonConvWrapper{F},
  variables::AbstractVector{<:DFGVariable},
  sfidx::Int,
  N::Integer;
  measurement = Vector{Tuple{}}(),
  needFreshMeasurements::Bool = true,
  solveKey::Symbol = :default,
) where {F <: AbstractFactor} # F might be Mixture
  #
  if length(variables) !== 0
    nothing
  else
    @debug("cannot prep ccw.param list with length(variables)==0, see DFG #590")
  end

  # in forward solve case, important to set which variable is being solved early in this sequence
  # set the 'solvefor' variable index -- i.e. which connected variable of the factor is being computed in this convolution. 
  ccwl.varidx[] = sfidx
  # ccwl.varidx[] = findfirst(==(solvefor), getLabel.(variables))
  
  # splice, type stable
  # make deepcopy of destination variable since multiple approxConv type computations should happen from different factors to the same variable
  tvarv = tuple(
    map(s->getVal(s; solveKey), variables[1:ccwl.varidx[]-1])..., 
    deepcopy(getVal(variables[ccwl.varidx[]]; solveKey)), # deepcopy(ccwl.varValsAll[][sfidx]),
    map(s->getVal(s; solveKey), variables[ccwl.varidx[]+1:end])...,
  )
  ccwl.varValsAll[] = tvarv
  
  # TODO, maxlen should parrot N (barring multi-/nullhypo issues)
  # everybody use maxlen number of points in belief function estimation
  maxlen = maximum((N, length.(ccwl.varValsAll[])...,))

  # if solving for more or less points in destination
  if N != length(ccwl.varValsAll[][ccwl.varidx[]])
    varT = getVariableType(variables[ccwl.varidx[]])
    # make vector right length
    resize!(ccwl.varValsAll[][ccwl.varidx[]], N)
    # define any new memory that might have been allocated
    for i in 1:N
      if !isdefined(ccwl.varValsAll[][ccwl.varidx[]], i)
        ccwl.varValsAll[][ccwl.varidx[]][i] = getPointDefault(varT)
      end
    end
  end

  # FIXME, confirm what happens when this is a partial dimension factor?  See #1246
  # indexing over all possible hypotheses
  xDim = getDimension(getVariableType(variables[ccwl.varidx[]]))
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

  return maxlen
end

function _beforeSolveCCW!(
  F_::Type{<:AbstractPrior},
  ccwl::CommonConvWrapper{F},
  variables::AbstractVector{<:DFGVariable},
  sfidx::Int,
  N::Integer;
  measurement = Vector{Tuple{}}(),
  needFreshMeasurements::Bool = true,
  solveKey::Symbol = :default,
) where {F <: AbstractFactor} # F might be Mixture
  # FIXME, NEEDS TO BE CLEANED UP AND WORK ON MANIFOLDS PROPER

  ccwl.varidx[] = sfidx
  @assert ccwl.varidx[] == 1 "Solving on Prior with CCW should have sfidx=1, priors are unary factors."

  # setup the partial or complete decision variable dimensions for this ccwl object
  # NOTE perhaps deconv has changed the decision variable list, so placed here during consolidation phase
  _setCCWDecisionDimsConv!(ccwl, getDimension(getVariableType(variables[ccwl.varidx[]])))

  solveForPts = getVal(variables[ccwl.varidx[]]; solveKey)
  maxlen = maximum([N; length(solveForPts); length(ccwl.varValsAll[][ccwl.varidx[]])])  # calcZDim(ccwl); length(measurement[1])

  # FIXME do not divert Mixture for sampling
  # update ccwl.measurement values
  updateMeasurement!(ccwl, maxlen; needFreshMeasurements, measurement, _allowThreads=true)

  return maxlen
end

# TODO, can likely deprecate this
function _beforeSolveCCW!(
  ccwl::Union{CommonConvWrapper{F}, CommonConvWrapper{Mixture{N_, F, S, T}}},
  Xi::AbstractVector{<:DFGVariable},
  # destVarVals::AbstractVector,
  sfidx::Int,
  N::Integer;
  kw...,
) where {N_, F <: AbstractRelative, S, T}
  #
  return _beforeSolveCCW!(F, ccwl, Xi, sfidx, N; kw...)
end

function _beforeSolveCCW!(
  ccwl::Union{CommonConvWrapper{F}, CommonConvWrapper{Mixture{N_, F, S, T}}},
  Xi::AbstractVector{<:DFGVariable},
  # destVarVals::AbstractVector,
  sfidx::Int,
  N::Integer;
  kw...,
) where {N_, F <: AbstractPrior, S, T}
  #
  return _beforeSolveCCW!(F, ccwl, Xi, sfidx, N; kw...)
end


#
