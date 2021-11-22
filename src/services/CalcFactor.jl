# New factor interface, something perhaps like this


export calcFactorResidualTemporary


getFactorOperationalMemoryType(dfg::SolverParams) = CommonConvWrapper
getFactorOperationalMemoryType(dfg::NoSolverParams) = CommonConvWrapper




function _getDimensionsPartial(ccw::CommonConvWrapper)
  # @warn "_getDimensionsPartial not ready for use yet"
  ccw.partialDims
end
_getDimensionsPartial(data::GenericFunctionNodeData) = _getCCW(data) |> _getDimensionsPartial
_getDimensionsPartial(fct::DFGFactor) = _getDimensionsPartial(_getCCW(fct))
_getDimensionsPartial(fg::AbstractDFG, lbl::Symbol) = _getDimensionsPartial(getFactor(fg, lbl))


"""
    $SIGNATURES
Get `.factormetadata` for each CPT in CCW for a specific factor in `fg`. 
"""
_getFMdThread(ccw::CommonConvWrapper, 
              thrid::Int=Threads.threadid()) = ccw.cpt[thrid].factormetadata
#
_getFMdThread(fc::Union{GenericFunctionNodeData,DFGFactor}, 
              thrid::Int=Threads.threadid()) = _getFMdThread(_getCCW(fc), thrid)
#
_getFMdThread(dfg::AbstractDFG,
              lbl::Symbol,
              thrid::Int=Threads.threadid()) = _getFMdThread(_getCCW(dfg, lbl), thrid)
#




# Helper function to construct CF from a CCW
function CalcFactor(ccwl::CommonConvWrapper;
                    factor = ccwl.usrfnc!, 
                    metadata = _getFMdThread(ccwl), 
                    _sampleIdx = 0, 
                    _measCount = length(ccwl.measurement),
                    _legacyMeas = ccwl.measurement, 
                    _legacyParams = ccwl.params, 
                    _allowThreads = true )
  #
  CalcFactor( factor,
              metadata,
              _sampleIdx,
              _measCount,
              _legacyMeas,
              _legacyParams,
              _allowThreads )
end


"""
    $SIGNATURES

Sample the factor stochastic model `N::Int` times and store the samples in the preallocated `ccw.measurement` container.

DevNotes
- Use in place operations where possible and remember `measurement` is a `::Tuple`.
- TODO only works on `.threadid()==1` at present, see #1094
- Also see, JuliaRobotics/RoME.jl#465
"""
sampleFactor(cf::CalcFactor{<:AbstractFactor}, N::Int=1  ) = [getSample(cf) for _=1:N]



function Base.show(io::IO, x::CalcFactor)
  println(io, )
  printstyled(io, " CalcFactor:\n", color=:blue)
  println(io, "  .factor: ", typeof(x.factor))
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
  try
    M = getManifold(T)
    return manifold_dimension(M)
  catch
    try 
      M = getManifold(cf.factor)
      return manifold_dimension(M)
    catch
      @warn "no method getManifold(::$(string(T))), calcZDim will attempt legacy length(sample) method instead"
    end
  end
  
  # NOTE try to make sure we get matrix back (not a vector)
  smpls = sampleFactor(cf, 2)[1]
  return length(smpls)
end

calcZDim(ccw::CommonConvWrapper) = calcZDim(CalcFactor(ccw))

calcZDim(cf::CalcFactor{<:GenericMarginal}) = 0

calcZDim(cf::CalcFactor{<:ManifoldPrior}) = manifold_dimension(cf.factor.M)



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
calcFactorResidual(dfgfct::DFGFactor, args...; ccw::CommonConvWrapper=IIF._getCCW(dfgfct)) = CalcFactor(ccw)(args...)
calcFactorResidual(dfg::AbstractDFG, fctsym::Symbol, args...) = calcFactorResidual(getFactor(dfg, fctsym), args...)


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
function calcFactorResidualTemporary( fct::AbstractRelative, 
                                      varTypes::Tuple,
                                      measurement,
                                      pts::Tuple;
                                      tfg::AbstractDFG = initfg(),
                                      _blockRecursion::Bool=false,
                                      doTime::Bool=false  )
  #

  # build a new temporary graph
  _, _dfgfct = _buildGraphByFactorAndTypes!(fct, varTypes, pts, dfg=tfg, _blockRecursion=_blockRecursion)
  
  # get a fresh measurement if needed
  _measurement = if measurement != [] #length(measurement) != 0
    measurement
  else
    # now use the CommonConvWrapper object in `_dfgfct`
    cfo = CalcFactor(_getCCW(_dfgfct))
    sampleFactor(cfo, 1)[1]
  end
  
  # assume a single sample point is being run
  return if doTime
    @time res = calcFactorResidual(_dfgfct, _measurement, pts...)
    res
  else
    calcFactorResidual(_dfgfct, _measurement, pts...)
  end
end


## =============================================================================================
## FactorOperationalMemory helper constructors
## =============================================================================================


function ConvPerThread( X::AbstractVector{P},
                        zDim::Int,
                        factormetadata::FactorMetadata;
                        particleidx::Int=1,
                        activehypo= 1:length(params),
                        perturb=zeros(zDim),
                        res=zeros(zDim),
                        thrid_ = 0  ) where P
  #
  return ConvPerThread{typeof(res), typeof(factormetadata), Any}( thrid_,
                        particleidx,
                        factormetadata,
                        Int[activehypo;],
                        perturb,
                        X,
                        res )
end



function CommonConvWrapper( fnc::T,
                            X::AbstractVector{P},
                            zDim::Int,
                            params::AbstractVector{<:AbstractVector{Q}},
                            factormetadata::FactorMetadata;
                            partial::Bool=false,
                            hypotheses::H=nothing,
                            certainhypo=nothing,
                            activehypo= 1:length(params),
                            nullhypo::Real=0,
                            varidx::Int=1,
                            measurement::AbstractVector=Vector(Vector{Float64}(),),
                            particleidx::Int=1,
                            xDim::Int=size(X,1),
                            partialDims::AbstractVector{<:Integer}=1:length(X),
                            perturb=zeros(zDim),
                            res::AbstractVector{<:Real}=zeros(zDim),
                            threadmodel::Type{<:_AbstractThreadModel}=MultiThreaded,
                            inflation::Real=3.0,
                            vartypes=typeof.(getVariableType.(factormetadata.fullvariables)),
                            gradients=nothing) where {T<:AbstractFactor,P,H,Q}
  #
  # FIXME, deprecate params::Vector throughout codebase
  tup = tuple(params...)
  nms = tuple(factormetadata.variablelist...)
  p_ntup = NamedTuple{nms,typeof(tup)}(tup)

  return  CommonConvWrapper(fnc,
                            xDim,
                            zDim,
                            partial,
                            hypotheses,
                            certainhypo,
                            Float64(nullhypo),
                            p_ntup, # params,
                            varidx,
                            measurement,
                            threadmodel,
                            (i->ConvPerThread(X, zDim,factormetadata, particleidx=particleidx,
                                              activehypo=activehypo, 
                                              perturb=perturb, res=res )).(1:Threads.nthreads()),
                            inflation,
                            partialDims,
                            DataType[vartypes...],
                            gradients)
end


function _resizePointsVector!(vecP::AbstractVector{P}, mkd::ManifoldKernelDensity, N::Int) where P
  #
  pN = length(vecP)
  resize!(vecP, N)
  for j in pN:N
    smp = AMP.sample(mkd, 1)[1]
    # @show j, smp, typeof(smp), typeof(vecP[j])
    vecP[j] = smp[1]
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

DevNotes
- FIXME ARR internally should become a NamedTuple
"""
function prepareparamsarray(Xi::Vector{<:DFGVariable},
                            solvefor::Union{Nothing, Symbol},
                            N::Int=0;
                            solveKey::Symbol=:default  ) where P
  #
  # FIXME ON FIRE, refactor to new NamedTuple instead
  varParamsAll = Vector{Vector{Any}}()

  LEN = Int[]
  maxlen = N # FIXME see #105
  count = 0
  sfidx = 0

  for xi in Xi
    vecP = getVal(xi, solveKey=solveKey)
    push!(varParamsAll, vecP)
    LEN = length.(varParamsAll)
    maxlen = maximum([N; LEN])
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
      _resizePointsVector!(varParamsAll[i], Pr, maxlen)
    end
  end

  # TODO --rather define reusable memory for the proposal
  # we are generating a proposal distribution, not direct replacement for existing memory and hence the deepcopy.
  if sfidx > 0 
    varParamsAll[sfidx] = deepcopy(varParamsAll[sfidx]) 
  end

  # get solvefor manifolds
  # FIXME deprecate use of (:null,)
  mani = length(Xi)==0 || sfidx==0 ? (:null,) : getManifold(Xi[sfidx])

  # FIXME, forcing maxlen to N results in errors (see test/testVariousNSolveSize.jl) see #105
  # maxlen = N == 0 ? maxlen : N
  return varParamsAll, maxlen, sfidx, mani
end

"""
    $SIGNATURES
Internal method to set which dimensions should be used as the decision variables for later numerical optimization.
"""
function _setCCWDecisionDimsConv!(ccwl::Union{CommonConvWrapper{F},
                                              CommonConvWrapper{Mixture{N_,F,S,T}}} ) where {N_,F<:Union{AbstractManifoldMinimize, AbstractRelativeMinimize, AbstractRelativeRoots, AbstractPrior},S,T}
  #

  # NOTE should only be done in the constructor
  ccwl.partialDims = if ccwl.partial
    Int[ccwl.usrfnc!.partial...]
  else
    Int[1:ccwl.xDim...]
  end
  
  nothing
end



"""
    $SIGNATURES

Notes
- Can be called with `length(Xi)==0`
"""
function _prepCCW(Xi::Vector{<:DFGVariable},
                  usrfnc::T;
                  multihypo::Union{Nothing, Distributions.Categorical}=nothing,
                  nullhypo::Real=0.0,
                  threadmodel=MultiThreaded,
                  inflation::Real=0.0,
                  _blockRecursion::Bool=false  ) where {T <: AbstractFactor}
  #
  length(Xi) !== 0 ? nothing : @debug("cannot prep ccw.param list with length(Xi)==0, see DFG #590")

  pttypes = getVariableType.(Xi) .|> getPointType
  PointType = 0 < length(pttypes) ? pttypes[1] : Vector{Float64}
  # FIXME stop using Any, see #1321
  # varParamsAll = Vector{Vector{Any}}()
  varParamsAll, maxlen, sfidx, mani = prepareparamsarray( Xi, nothing, 0) # Nothing for init.

  # standard factor metadata
  sflbl = 0==length(Xi) ? :null : getLabel(Xi[end])
  lbs = getLabel.(Xi)
  fmd = FactorMetadata(Xi, lbs, varParamsAll, sflbl, nothing)

  # create a temporary CalcFactor object for extracting the first sample
  # TODO, deprecate this:  guess measurement points type
  # MeasType = Vector{Float64} # FIXME use `usrfnc` to get this information instead
  _cf = CalcFactor( usrfnc, fmd, 0, 1, nothing, varParamsAll, false) # (Vector{MeasType}(),)
  
  # get a measurement sample
  meas_single = sampleFactor(_cf, 1)[1]

  #TODO preallocate measurement?
  measurement = Vector{typeof(meas_single)}()

  # get the measurement dimension
  zdim = calcZDim(_cf)
  # some hypo resolution
  certainhypo = multihypo !== nothing ? collect(1:length(multihypo.p))[multihypo.p .== 0.0] : collect(1:length(Xi))
  
  # partialDims are sensitive to both which solvefor variable index and whether the factor is partial
  ispartl = hasfield(T, :partial)
  partialDims = if ispartl
    Int[usrfnc.partial...]
  else
    Int[]
  end

  # as per struct CommonConvWrapper
  varTypes::Vector{DataType} = typeof.(getVariableType.(Xi))
  gradients = nothing
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
      pts_ = map(x->x[1], varParamsAll)
      # FIXME, only using first meas and params values at this time...
      # NOTE, must block recurions here, since FGC uses this function to calculate numerical gradients on a temp fg.
      # assume for now fractional-var in multihypo have same varType
      hypoidxs = _selectHypoVariables(pts_, multihypo)
      gradients = FactorGradientsCached!(usrfnc, tuple(varTypes[hypoidxs]...), measurement_, tuple(pts_[hypoidxs]...), _blockRecursion=true);
    end
  catch e
    @warn "Unable to create measurements and gradients for $usrfnc during prep of CCW, falling back on no-partial information assumption.  Enable ENV[\"JULIA_DEBUG\"] = \"IncrementalInference\" for @debug printing to see the error."
    # rethrow(e)
    @debug(e)
  end

  ccw = CommonConvWrapper(
          usrfnc,
          PointType[],
          zdim,
          varParamsAll,
          fmd;
          partial = ispartl,
          measurement,
          hypotheses=multihypo,
          certainhypo=certainhypo,
          nullhypo=nullhypo,
          threadmodel=threadmodel,
          inflation=inflation,
          partialDims=partialDims,
          vartypes = varTypes,
          gradients=gradients
        )
  #
  return ccw
end

"""
    $(SIGNATURES)

Prepare a common functor computation object `prepareCommonConvWrapper{T}` containing 
the user factor functor along with additional variables and information using during 
approximate convolution computations.

DevNotes
- TODO consolidate with others, see https://github.com/JuliaRobotics/IncrementalInference.jl/projects/6
"""
function prepareCommonConvWrapper!( F_::Type{<:AbstractRelative},
                                    ccwl::CommonConvWrapper{F},
                                    Xi::AbstractVector{<:DFGVariable},
                                    solvefor::Symbol,
                                    N::Int;
                                    needFreshMeasurements::Bool=true,
                                    solveKey::Symbol=:default  ) where {F <: AbstractFactor}
  #

  # FIXME, order of fmd ccwl cf are a little weird and should be revised.
  pttypes = getVariableType.(Xi) .|> getPointType
  PointType = 0 < length(pttypes) ? pttypes[1] : Vector{Float64}

  #FIXME, see #1321
  # vecPtsArr = Vector{Vector{Any}}()

  #TODO some better consolidate is needed
  ccwl.vartypes = typeof.(getVariableType.(Xi))

  # FIXME maxlen should parrot N (barring multi-/nullhypo issues)
  vecPtsArr, maxlen, sfidx, mani = prepareparamsarray( Xi, solvefor, N, solveKey=solveKey)

  # FIXME ON FIRE, what happens if this is a partial dimension factor?  See #1246
  ccwl.xDim = getDimension(getVariableType(Xi[sfidx]))
  # ccwl.xDim = length(vecPtsArr[sfidx][1])
  # TODO should be selecting for the correct multihypothesis mode

    # setup the partial or complete decision variable dimensions for this ccwl object
    # NOTE perhaps deconv has changed the decision variable list, so placed here during consolidation phase
    # TODO, should this not be part of `prepareCommonConvWrapper` -- only here do we look for .partial
    _setCCWDecisionDimsConv!(ccwl)

  # SHOULD WE SLICE ARR DOWN BY PARTIAL DIMS HERE (OR LATER)?
  # FIXME refactor new type higher up.
  tup = tuple(vecPtsArr...)
  nms = tuple(getLabel.(Xi)...)
  ntp = NamedTuple{nms,typeof(tup)}(tup)
  ccwl.params = ntp
  
  # get factor metadata -- TODO, populate, also see #784
  fmd = FactorMetadata(Xi, getLabel.(Xi), ccwl.params, solvefor, nothing)

  # TODO consolidate with ccwl??
  # FIXME do not divert Mixture for sampling
  # cf = _buildCalcFactorMixture(ccwl, fmd, 1, ccwl.measurement, ccwl.params) # TODO perhaps 0 is safer
  cf = CalcFactor( ccwl.usrfnc!, fmd, 0, length(ccwl.measurement), ccwl.measurement, ccwl.params, true)

  #  get variable node data
  vnds = Xi

  # option to disable fresh samples
  if needFreshMeasurements
    # TODO refactor
    ccwl.measurement = sampleFactor(cf, maxlen)
    # sampleFactor!(ccwl, maxlen, fmd, vnds)
  end

  # cache the measurement dimension
  ccwl.zDim = calcZDim(CalcFactor(ccwl))
  # set the 'solvefor' variable index -- i.e. which connected variable of the factor is being computed in this convolution. 
  ccwl.varidx = sfidx

  # set each CPT
  for thrid in 1:Threads.nthreads()
    cpt_ = ccwl.cpt[thrid] 
    cpt_.X = ccwl.params[sfidx]

    # used in ccw functor for AbstractRelativeMinimize
    # TODO JT - Confirm it should be updated here. Testing in _prepCCW
    resize!(cpt_.res, ccwl.zDim) 
    fill!(cpt_.res, 0.0)
  end

  # calculate new gradients perhaps
  # J = ccwl.gradients(measurement..., pts...)

  return sfidx, maxlen, mani
end


function prepareCommonConvWrapper!( ccwl::Union{CommonConvWrapper{F},
                                                CommonConvWrapper{Mixture{N_,F,S,T}}},
                                    Xi::AbstractVector{<:DFGVariable},
                                    solvefor::Symbol,
                                    N::Int;
                                    kw...  ) where {N_,F<:AbstractRelative,S,T}
  #
  prepareCommonConvWrapper!(F, ccwl, Xi, solvefor, N; kw...)
end






#