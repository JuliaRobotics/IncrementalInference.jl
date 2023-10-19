
# TODO deprecate testshuffle
function _checkErrorCCWNumerics(
  ccwl::Union{CommonConvWrapper{F}, CommonConvWrapper{Mixture{N_, F, S, T}}},
  testshuffle::Bool = false,
) where {N_, F <: AbstractRelativeMinimize, S, T}
  return nothing
end
function _checkErrorCCWNumerics(
  ccwl::Union{CommonConvWrapper{F}, CommonConvWrapper{Mixture{N_, F, S, T}}},
  testshuffle::Bool = false,
) where {N_, F <: AbstractManifoldMinimize, S, T}
  return nothing
end


function _perturbIfNecessary(
  fcttype::Union{F, <:Mixture{N_, F, S, T}},
  len::Int = 1,
  perturbation::Real = 1e-10,
) where {N_, F <: AbstractRelativeMinimize, S, T}
  return 0
end

function _perturbIfNecessary(
  fcttype::Union{F, <:Mixture{N_, F, S, T}},
  len::Int = 1,
  perturbation::Real = 1e-10,
) where {N_, F <: AbstractManifoldMinimize, S, T}
  return 0
end
#


# internal use only, and selected out from approxDeconv functions
function _solveLambdaNumeric(
  fcttype::AbstractPrior,
  objResX::Function,
  residual::AbstractVector{<:Real},
  u0::AbstractVector{<:Real},
  islen1::Bool = false;
  perturb::Real = 1e-10,
)
  return u0
end
#


function _solveLambdaNumeric(
  fcttype::Union{F, <:Mixture{N_, F, S, T}},
  objResX::Function,
  residual::AbstractVector{<:Real},
  u0::AbstractVector{<:Real},
  islen1::Bool = false,
) where {N_, F <: AbstractRelativeMinimize, S, T}
  # retries::Int=3 )
  #
  # wrt #467 allow residual to be standardize for Roots and Minimize and Parametric cases.
  r = if islen1
    Optim.optimize((x) -> (residual .= objResX(x); sum(residual .^ 2)), u0, Optim.BFGS())
  else
    Optim.optimize((x) -> (residual .= objResX(x); sum(residual .^ 2)), u0, Optim.Options(;iterations=1000))
  end


  if !Optim.converged(r)
    @warn "Optim did not converge:" r maxlog=10
  end

  # 
  return r.minimizer
end

struct OptimCalcConv end
# CalcFactorNormSq cost function for an input in coordinates as used by Optim.jl
function (hypoCalcFactor::CalcFactorNormSq)(::Type{OptimCalcConv}, M::AbstractManifold, Xc::AbstractVector)
  # hypoCalcFactor.manifold is the factor's manifold, not the variable's manifold that is needed here
  ϵ = getPointIdentity(M)
  X = get_vector(M, ϵ, SVector(Xc), DefaultOrthogonalBasis())
  p = exp(M, ϵ, X)
  return hypoCalcFactor(CalcConv, p)
end
(hypoCalcFactor::CalcFactorNormSq)(M::AbstractManifold, p) = hypoCalcFactor(OptimCalcConv, M, p)

struct ManoptCalcConv end

function (hypoCalcFactor::CalcFactorNormSq)(::Type{ManoptCalcConv}, M::AbstractManifold, p)
  return hypoCalcFactor(CalcConv, p)
end

function _solveLambdaNumeric(
  fcttype::Union{F, <:Mixture{N_, F, S, T}},
  hypoCalcFactor,
  residual::AbstractVector{<:Real},
  u0,#::AbstractVector{<:Real},
  variableType::InferenceVariable,
  islen1::Bool = false,
) where {N_, F <: AbstractManifoldMinimize, S, T}
  #
  M = getManifold(variableType)
  # the variable is a manifold point, we are working on the tangent plane in optim for now.
  # 
  #TODO this is not general to all manifolds, should work for lie groups.
  ϵ = getPointIdentity(variableType)

  X0c = zero(MVector{getDimension(M),Float64})
  X0c .= vee(M, u0, log(M, ϵ, u0))

  alg = islen1 ? Optim.BFGS() : Optim.NelderMead()

  #WIP extremely slow, but runs, mean in manopt is bottleneck
  # just to show how we can now swop to manopt
  if false
    r = Manopt.NelderMead(
      M, 
      (M, x)->hypoCalcFactor(ManoptCalcConv, M, x),
      NelderMeadSimplex(M, u0, DefaultOrthogonalBasis());
      retraction_method = ExponentialRetraction()
    )
    return r
  elseif false
    r = gradient_descent(
      M,
      (M,x)->hypoCalcFactor(x),
      (M, x)-> factorGradient(hypoCalcFactor, M, x),
      u0;
      stepsize=ConstantStepsize(0.1), 
    )
    return r
  end

  r = Optim.optimize(
    x->hypoCalcFactor(OptimCalcConv, M, x),
    X0c,
    alg
  )
  
  if !Optim.converged(r)
    # TODO find good way for a solve to store diagnostics about number of failed converges etc.
    @warn "Optim did not converge (maxlog=10):" r maxlog=10
  end
  return exp(M, ϵ, hat(M, ϵ, r.minimizer))
end

## deconvolution with calcfactor wip
struct CalcDeconv end

function (cf::CalcFactorNormSq)(::Type{CalcDeconv}, meas) 
  res = cf(meas, map(vvh -> _getindex_anyn(vvh, cf._sampleIdx), cf._legacyParams)...)
  return sum(x->x^2, res)
end

# for deconv with the measurement a tangent vector, can dispatch for other measurement types.
function (hypoCalcFactor::CalcFactorNormSq)(::Type{CalcDeconv}, M::AbstractManifold, Xc::AbstractVector)
  ϵ = getPointIdentity(M)
  X = get_vector(M, ϵ, Xc, DefaultOrthogonalBasis())
  return hypoCalcFactor(CalcDeconv, X)
end

# NOTE Optim.jl version that assumes measurement is on the tangent
# TODO test / dev for n-ary factor deconv
# TODO Consolidate with _solveLambdaNumeric, see #1374
function _solveLambdaNumericMeas(
  fcttype::Union{F, <:Mixture{N_, F, S, T}},
  hypoCalcFactor,
  X0,#::AbstractVector{<:Real},
  islen1::Bool = false,
) where {N_, F <: AbstractManifoldMinimize, S, T}
  #
  M = getManifold(fcttype)
  ϵ = getPointIdentity(M)
  X0c = zeros(manifold_dimension(M))
  X0c .= vee(M, ϵ, X0)

  alg = islen1 ? Optim.BFGS() : Optim.NelderMead()

  r = Optim.optimize(
    x->hypoCalcFactor(CalcDeconv, M, x),
    X0c,
    alg
  )
  if !Optim.converged(r)
    @debug "Optim did not converge:" r
  end

  return hat(M, ϵ, r.minimizer)
end

## ================================================================================================
## Heavy dispatch for all AbstractFactor / Mixture cases below
## ================================================================================================

# internal function to dispatch view on either vector or matrix, rows are dims and samples are columns
_getindextuple(tup::Tuple, ind1::Int) = [getindex(t, ind1) for t in tup]

_getusrfnc(ccwl::CommonConvWrapper) = ccwl.usrfnc!
_getusrfnc(ccwl::CommonConvWrapper{<:Mixture}) = ccwl.usrfnc!.mechanics

function _buildCalcFactor(
  ccwl::CommonConvWrapper,
  smpid,
  varParams,
  activehypo,
  _slack = nothing,
)
  #
  # FIXME, make thread safe (cache)
  # activevariables = view(ccwl.fullvariables, activehypo)
  activevariables = ccwl.fullvariables[activehypo]

  solveforidx = findfirst(==(ccwl.varidx[]), activehypo)

  return CalcFactorNormSq(
    _getusrfnc(ccwl),           #factor
    smpid,                      #_sampleIdx
    varParams,                  #_legacyParams
    true,                       #_allowThreads
    ccwl.dummyCache,            #_cache
    tuple(activevariables...),  #fullvariables
    solveforidx,                #solvefor
    getManifold(ccwl),           #manifold
    ccwl.measurement,
    _slack,
  )
end

"""
    $SIGNATURES
Internal function to build lambda pre-objective function for finding factor residuals. 

DevNotes
- TODO refactor relationship and common fields between (CCW, FMd, CPT, CalcFactor)
"""
function _buildCalcFactorLambdaSample(
  ccwl::CommonConvWrapper,
  smpid::Integer,
  measurement_, # since JL@v1.9, don't use default ccwl.measurement here, must pass from caller
  _slack = nothing,
)
  #

  # TODO from obsolete _view:
  # Should be replaced with ccw.hypoParams::Tuple(hypo1, hypo2,...), made at construction and allows direct hypo lookup
  # DevNotes, also see new `hyporecipe` approach (towards consolidation CCW CPT FMd CF...)

  # build a view to the decision variable memory
  varValsHypo = ccwl.varValsAll[][ccwl.hyporecipe.activehypo]

  # get the operational CalcFactor object
  cf = _buildCalcFactor(ccwl, smpid, varValsHypo, ccwl.hyporecipe.activehypo)

  # reset the residual vector
  fill!(ccwl.res, 0.0) # Roots->xDim | Minimize->zDim

  # build static lambda
  unrollHypo! = if _slack === nothing
    # DESIGN DECISION WAS MADE THAT CALCFACTOR CALLS DO NOT DO INPLACE CHANGES TO ARGUMENTS, INSTEAD USING ISBITSTYPEs!!!!!!!!!
    # 5.366727 seconds (17.48 M allocations: 893.768 MiB, 8.76% gc time)
    # () -> (cf::CalcFactorNormSq)(measurement_, smpid, varValsHypo)
    # 6.075632 seconds (19.73 M allocations: 919.118 MiB, 9.14% gc time)
    () -> cf(measurement_[smpid], map(vvh -> _getindex_anyn(vvh, smpid), varValsHypo)...)
  else
    # slack is used to shift the residual away from the natural "zero" tension position of a factor, 
    # this is useful when calculating factor gradients at a variety of param locations resulting in "non-zero slack" of the residual.
    # see `IIF.calcFactorResidualTemporary`
    # NOTE this minus operation assumes _slack is either coordinate or tangent vector element (not a manifold or group element)
    () ->
      cf(measurement_[smpid], map(vvh -> _getindex_anyn(vvh, smpid), varValsHypo)...) .- _slack
  end

  return unrollHypo!
end

"""
    $(SIGNATURES)

Solve free variable x by root finding residual function `fgr.usrfnc(res, x)`.  This is the 
penultimate step before calling numerical operations to move actual estimates, which is 
done by an internally created lambda function.

Notes
- Assumes `cpt_.p` is already set to desired X decision variable dimensions and size. 
- Assumes only `ccw.particleidx` will be solved for
- small random (off-manifold) perturbation used to prevent trivial solver cases, div by 0 etc.
  - perturb is necessary for NLsolve (obsolete) cases, and smaller than 1e-10 will result in test failure
- Also incorporates the active hypo lookup

DevNotes
- TODO testshuffle is now obsolete, should be removed
- TODO perhaps consolidate perturbation with inflation or nullhypo
"""
function _solveCCWNumeric!(
  ccwl::Union{<:CommonConvWrapper{F}, <:CommonConvWrapper{<:Mixture{N_, F, S, T}}},
  _slack = nothing;
  perturb::Real = 1e-10,
) where {N_, F <: AbstractRelative, S, T}
  #
  
  #
  # thrid = Threads.threadid()
  smpid = ccwl.particleidx[]
  # cannot Nelder-Mead on 1dim, partial can be 1dim or more but being conservative.
  islen1 = length(ccwl.partialDims) == 1 || ccwl.partial
  # islen1 = length(cpt_.X[:, smpid]) == 1 || ccwl.partial

  # NOTE the factor residual function will receive as input args a slice from ccwl.varValsAll, hence 
  #  ccwl.varValsAll[][ccwl.varidx[]] and target should point to the same memory; BUT
  #  remember that during approxConv the graph variable cannot be directly updated and
  #  a separate deepcopy of the destination (aka target) memory is necessary.
  #  Choosen solution is to splice together ccwl.varValsAll each time, with destination as 
  #  deepcopy but other input variables are just point to the source variable values directly. 
  target = if ccwl.partial  # FIXME likely type-instability on `typeof(target)`
    # view(ccwl.varValsAll[][ccwl.varidx[]][smpid], ccwl.partialDims)
    ccwl.varValsAll[][ccwl.varidx[]][smpid][ccwl.partialDims]
  else
    ccwl.varValsAll[][ccwl.varidx[]][smpid]
  end
  # build the pre-objective function for this sample's hypothesis selection
  unrollHypo! = _buildCalcFactorLambdaSample(
      # destVarVals,
      ccwl, 
      smpid,
      ccwl.measurement,
      _slack,
  )

  # broadcast updates original view memory location
  ## using CalcFactor legacy path inside (::CalcFactor)

  # _hypoObj = (x) -> (target[] = x; unrollHypo!())
  function _hypoObj(x)
    copyto!(target, x)
    return unrollHypo!()
  end

  # TODO small off-manifold perturbation is a numerical workaround only, make on-manifold requires RoME.jl #244
  # use all element dimensions : ==> 1:ccwl.xDim
  # target .+= _perturbIfNecessary(getFactorType(ccwl), length(target), perturb)
  sfidx = ccwl.varidx[]
  # do the parameter search over defined decision variables using Minimization
  X = ccwl.varValsAll[][sfidx][smpid][ccwl.partialDims]
  # X = if ccwl.partial # TODO check for type-instability on `X`
  #   collect(view(ccwl.varValsAll[][sfidx][smpid], ccwl.partialDims))
  # else
  #   ccwl.varValsAll[][sfidx][smpid][ccwl.partialDims]
  # end
  # # X = destVarVals[smpid]#[ccwl.partialDims]
      
  retval = _solveLambdaNumeric(
    getFactorType(ccwl), 
    _hypoObj, 
    ccwl.res, 
    X, 
    islen1
  )

  # Check for NaNs
  if sum(isnan.(retval)) != 0
    @error "$(ccwl.usrfnc!), got NaN, smpid = $(smpid), r=$(retval)\n"
    return nothing
  end

  # insert result back at the correct variable element location
  if ccwl.partial
    # NOTE use workaround of TranslationGroup for coordinates on partial assignment
    # FIXME consolidate to Manopt and upgrade to Riemannian (i.e. incl non-groups)
    M = getManifold(ccwl) # TranslationGroup(length(ccwl.varValsAll[][sfidx][smpid]))
    src = Vector{typeof(retval)}()
    push!(src, retval)
    setPointPartial!(M, ccwl.varValsAll[][sfidx], M, src, ccwl.partialDims, smpid, 1, true )
    # ccwl.varValsAll[][sfidx][smpid][ccwl.partialDims] .= retval
  else
    # copyto!(ccwl.varValsAll[sfidx][smpid], retval)
    copyto!(ccwl.varValsAll[][sfidx][smpid][ccwl.partialDims], retval)
  end

  return nothing
end
# brainstorming
# should only be calling a new arg list according to activehypo at start of particle
# Try calling an existing lambda
# sensitive to which hypo of course , see #1024
#

struct CalcConv end

_getindex_anyn(vec, n) = begin
  len = length(vec)
  # 1:len or any random element in that range
  getindex(vec, n <= len ? n : rand(1:len) )
end

# NOTE to future self, this will likely become the cost function for Manopt as:
# function (cf::CalcFactorNormSq)(M::AbstractManifold, x)
# CalcConv is likeley needed for conv vs deconv
function (cf::CalcFactorNormSq)(::Type{CalcConv}, x) 
  
  sampleIdx = cf._sampleIdx
  varValsHypo = cf._legacyParams
  # set the target hypo on the correct sample to free variable x, was target object
  varValsHypo[cf.solvefor][sampleIdx] = x

  res = cf(cf.measurement[sampleIdx], map(vvh -> _getindex_anyn(vvh, sampleIdx), varValsHypo)...)
  res = isnothing(cf.slack) ? res : res .- cf.slack
  return sum(x->x^2, res)
end
#default to conv
(cf::CalcFactorNormSq)(x) = cf(CalcConv, x)

function _buildHypoCalcFactor(ccwl::CommonConvWrapper, smpid::Integer, _slack=nothing)
  # build a view to the decision variable memory
  varValsHypo = ccwl.varValsAll[][ccwl.hyporecipe.activehypo]
  # create calc factor selected hypo and samples
  #TODO lots of allocations, can we refactor to reuse?
  cf = _buildCalcFactor(
    ccwl,                        #
    smpid,                       # ends in _sampleIdx
    varValsHypo,                 # ends in _legacyParams
    ccwl.hyporecipe.activehypo,   # ends in solvefor::Int
    _slack,
  )
  return cf
end

function _solveCCWNumeric!(
  ccwl::Union{<:CommonConvWrapper{F}, <:CommonConvWrapper{<:Mixture{N_, F, S, T}}},
  _slack = nothing;
  perturb::Real = 1e-10,
) where {N_, F <: AbstractManifoldMinimize, S, T}
  #
  #   # FIXME, move this check higher and out of smpid loop
  # _checkErrorCCWNumerics(ccwl, testshuffle)

  smpid = ccwl.particleidx[]
  # cannot Nelder-Mead on 1dim, partial can be 1dim or more but being conservative.
  islen1 = length(ccwl.partialDims) == 1 || ccwl.partial

  # build the pre-objective function for this sample's hypothesis selection

  # SUPER IMPORTANT, this `target` is mem pointer that will be updated by optim library
  # target = view(ccwl.varValsAll[][ccwl.varidx[]], smpid) 

  # SUPER IMPORTANT ON PARTIALS, RESIDUAL FUNCTION MUST DEAL WITH PARTIAL AND WILL GET FULL VARIABLE POINTS REGARDLESS
  _hypoCalcFactor = _buildHypoCalcFactor(ccwl, smpid, _slack)

  # do the parameter search over defined decision variables using Minimization
  sfidx = ccwl.varidx[]
  u0 = ccwl.varValsAll[][ccwl.varidx[]][smpid] # starting point for optimization
  retval = _solveLambdaNumeric(
    getFactorType(ccwl),
    _hypoCalcFactor,
    ccwl.res,
    u0,
    getVariableType(ccwl.fullvariables[sfidx]), # only used for getting variable manifold and identity_element
    islen1,
  )

  # TBD Check for NaNs

  # NOTE insert result back at the correct variable element location
  ccwl.varValsAll[][ccwl.varidx[]][smpid] = retval

  return nothing
end

#
