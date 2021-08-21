# New factor interface, something perhaps like this


export calcFactorResidualTemporary


getFactorOperationalMemoryType(dfg::SolverParams) = CommonConvWrapper
getFactorOperationalMemoryType(dfg::NoSolverParams) = CommonConvWrapper


# Helper function to construct CF from a CCW
CalcFactor(ccwl::CommonConvWrapper) = CalcFactor( ccwl.usrfnc!, _getFMdThread(ccwl), 0, length(ccwl.measurement), ccwl.measurement, ccwl.params)



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
  return length(smpls[1])
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
fg = generateCanonicalFG_Kaess()

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
                                      measurement::Vector{<:Tuple},
                                      pts::Tuple;
                                      tfg::AbstractDFG = initfg(),
                                      _blockRecursion::Bool=false )
  #

  # build a new temporary graph
  _, _dfgfct = _buildGraphByFactorAndTypes!(fct, varTypes, pts, dfg=tfg, _blockRecursion=_blockRecursion)
  
  # get a fresh measurement if needed
  _measurement = if length(measurement) != 0
    measurement
  else
    # now use the CommonConvWrapper object in `_dfgfct`
    cfo = CalcFactor(_getCCW(_dfgfct))
    sampleFactor(cfo, 1)
  end

  # assume a single sample point is being run
  return calcFactorResidual(_dfgfct, _measurement[1]..., pts...)
end


## =============================================================================================
## FactorOperationalMemory helper constructors
## =============================================================================================


function ConvPerThread( X::AbstractVector{P},
                        zDim::Int,
                        factormetadata::FactorMetadata;
                        particleidx::Int=1,
                        activehypo= 1:length(params),
                        p::AbstractVector{<:Integer}=collect(1:1),
                        perturb=zeros(zDim),
                        res=zeros(zDim),
                        thrid_ = 0  ) where P
  #
  return ConvPerThread{typeof(res), typeof(factormetadata), Any}( thrid_,
                        particleidx,
                        factormetadata,
                        Int[activehypo;],
                        Int[p...;],
                        perturb,
                        X,
                        res )
end



function CommonConvWrapper( fnc::T,
                            X::AbstractVector{P},
                            zDim::Int,
                            params::AbstractVector{<:AbstractVector{Q}},
                            factormetadata::FactorMetadata;
                            specialzDim::Bool=false,
                            partial::Bool=false,
                            hypotheses::H=nothing,
                            certainhypo=nothing,
                            activehypo= 1:length(params),
                            nullhypo::Real=0,
                            varidx::Int=1,
                            measurement::Vector{<:Tuple}=Vector(Vector{Float64}(),),  # FIXME should not be a Matrix
                            particleidx::Int=1,
                            xDim::Int=size(X,1),
                            partialDims::AbstractVector{<:Integer}=collect(1:size(X,1)), # TODO make this SVector, and name partialDims
                            perturb=zeros(zDim),
                            res::AbstractVector{<:Real}=zeros(zDim),
                            threadmodel::Type{<:_AbstractThreadModel}=MultiThreaded,
                            inflation::Real=3.0,
                            vartypes=typeof.(getVariableType.(factormetadata.fullvariables)),
                            gradients=nothing) where {T<:FunctorInferenceType,P,H,Q}
  #
  return  CommonConvWrapper(fnc,
                            xDim,
                            zDim,
                            specialzDim,
                            partial,
                            hypotheses,
                            certainhypo,
                            Float64(nullhypo),
                            params,
                            varidx,
                            measurement,
                            threadmodel,
                            (i->ConvPerThread(X, zDim,factormetadata, particleidx=particleidx,
                                              activehypo=activehypo, p=partialDims, 
                                              perturb=perturb, res=res )).(1:Threads.nthreads()),
                            inflation,
                            partialDims,  # SVector(Int32.()...)
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
      _resizePointsVector!(ARR[i], Pr, maxlen)
    end
  end

  # TODO --rather define reusable memory for the proposal
  # we are generating a proposal distribution, not direct replacement for existing memory and hence the deepcopy.
  if sfidx > 0 
    ARR[sfidx] = deepcopy(ARR[sfidx]) 
  end

  # get solvefor manifolds
  # FIXME deprecate use of (:null,)
  mani = length(Xi)==0 || sfidx==0 ? (:null,) : getManifold(Xi[sfidx])

  # FIXME, forcing maxlen to N results in errors (see test/testVariousNSolveSize.jl) see #105
  # maxlen = N == 0 ? maxlen : N
  return maxlen, sfidx, mani
end

"""
    $SIGNATURES
Internal method to set which dimensions should be used as the decision variables for later numerical optimization.
"""
function _setCCWDecisionDimsConv!(ccwl::Union{CommonConvWrapper{F},
                                              CommonConvWrapper{Mixture{N_,F,S,T}}} ) where {N_,F<:Union{AbstractManifoldMinimize, AbstractRelativeMinimize, AbstractRelativeRoots, AbstractPrior},S,T}
  #
  # return nothing

  p = if ccwl.partial
    Int32[ccwl.usrfnc!.partial...]
  else
    Int32[1:ccwl.xDim...]
  end

  ccwl.partialDims = (p)
  # NOTE should only be done in the constructor
  for thrid in 1:Threads.nthreads()
    length(ccwl.cpt[thrid].p) != length(p) ? resize!(ccwl.cpt[thrid].p, length(p)) : nothing
    ccwl.cpt[thrid].p .= p # SVector... , see ccw.partialDims
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
  vecPtsArr = Vector{Vector{Any}}()

  #TODO some better consolidate is needed
  ccwl.vartypes = typeof.(getVariableType.(Xi))

  # FIXME maxlen should parrot N (barring multi-/nullhypo issues)
  maxlen, sfidx, mani = prepareparamsarray!(vecPtsArr, Xi, solvefor, N, solveKey=solveKey)

  # FIXME ON FIRE, what happens if this is a partial dimension factor?  See #1246
  ccwl.xDim = getDimension(getVariableType(Xi[sfidx]))
  # ccwl.xDim = length(vecPtsArr[sfidx][1])
  # TODO should be selecting for the correct multihypothesis mode

    # setup the partial or complete decision variable dimensions for this ccwl object
    # NOTE perhaps deconv has changed the decision variable list, so placed here during consolidation phase
    # TODO, should this not be part of `prepareCommonConvWrapper` -- only here do we look for .partial
    _setCCWDecisionDimsConv!(ccwl)

  # SHOULD WE SLICE ARR DOWN BY PARTIAL DIMS HERE (OR LATER)?
  ccwl.params = vecPtsArr # map( ar->view(ar, ccwl.partialDims, :), vecPtsArr)
  
  # get factor metadata -- TODO, populate, also see #784
  fmd = FactorMetadata(Xi, getLabel.(Xi), ccwl.params, solvefor, nothing)

  # TODO consolidate with ccwl??
  # FIXME do not divert Mixture for sampling
  # cf = _buildCalcFactorMixture(ccwl, fmd, 1, ccwl.measurement, ccwl.params) # TODO perhaps 0 is safer
  cf = CalcFactor( ccwl.usrfnc!, fmd, 0, length(ccwl.measurement), ccwl.measurement, ccwl.params)

  #  get variable node data
  vnds = Xi

  # option to disable fresh samples
  if needFreshMeasurements
    # TODO refactor
    ccwl.measurement = sampleFactor(cf, maxlen)
    # sampleFactor!(ccwl, maxlen, fmd, vnds)
  end


  ccwl.zDim = calcZDim(CalcFactor(ccwl))
  # if ccwl.specialzDim
  #   ccwl.zDim = ccwl.usrfnc!.zDim[sfidx]
  # else
  # end
  ccwl.varidx = sfidx

  # set each CPT
  for thrid in 1:Threads.nthreads()
    cpt_ = ccwl.cpt[thrid] 
    cpt_.X = ccwl.params[sfidx]

    # used in ccw functor for AbstractRelativeMinimize
    # TODO JT - Confirm it should be updated here. Testing in prepgenericconvolution
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