
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

function _checkErrorCCWNumerics(
  ccwl::Union{CommonConvWrapper{F}, CommonConvWrapper{Mixture{N_, F, S, T}}},
  testshuffle::Bool = false,
) where {N_, F <: AbstractRelativeRoots, S, T}
  #
  # error("<:AbstractRelativeRoots is obsolete, use one of the other <:AbstractRelative types instead.")
  # TODO get xDim = getDimension(getVariableType(Xi[sfidx])) but without having Xi
  if testshuffle || ccwl.partial
    error(
      "<:AbstractRelativeRoots factors with less or more measurement dimensions than variable dimensions have been discontinued, rather use <:AbstractManifoldMinimize.",
    )
  # elseif !(_getZDim(ccwl) >= ccwl.xDim && !ccwl.partial)
  #   error("Unresolved numeric <:AbstractRelativeRoots solve case")
  end
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

function _perturbIfNecessary(
  fcttype::Union{F, <:Mixture{N_, F, S, T}},
  len::Int = 1,
  perturbation::Real = 1e-10,
) where {N_, F <: AbstractRelativeRoots, S, T}
  return perturbation * randn(len)
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
) where {N_, F <: AbstractRelativeRoots, S, T}
  #

  #
  r = NLsolve.nlsolve((res, x) -> res .= objResX(x), u0; inplace = true) #, ftol=1e-14)

  #
  return r.zero
end

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
    Optim.optimize((x) -> (residual .= objResX(x); sum(residual .^ 2)), u0)
  end

  # 
  return r.minimizer
end

function _solveLambdaNumeric(
  fcttype::Union{F, <:Mixture{N_, F, S, T}},
  objResX::Function,
  residual::AbstractVector{<:Real},
  u0,#::AbstractVector{<:Real},
  variableType::InferenceVariable,
  islen1::Bool = false,
) where {N_, F <: AbstractManifoldMinimize, S, T}
  #
  M = getManifold(variableType) #fcttype.M
  # the variable is a manifold point, we are working on the tangent plane in optim for now.
  # 
  #TODO this is not general to all manifolds, should work for lie groups.
  # ϵ = identity_element(M, u0)
  ϵ = getPointIdentity(variableType)
  # X0c = get_coordinates(M, u0, log(M, ϵ, u0), DefaultOrthogonalBasis()) 
  X0c = vee(M, u0, log(M, ϵ, u0))

  # objResX(p) returns tangent vector at p, X=log(M, p, ...)
  # norm(M, p, X) == distance(M, p, X)
  #TODO fix closure for performance
  fM = getManifold(fcttype)
  function cost(p, X, Xc)
    hat!(M, X, ϵ, Xc)
    retract!(M, p, ϵ, X)
    # X = objResX(p)
    # return norm(fM, p, X)^2 #TODO the manifold of p and X are not always the same
    #options getPointIdentity or leave it to factor 
    residual = objResX(p)
    return sum(residual .^ 2)
  end

  # # separate statements to try improve type-stability 
  # r = if islen1
  #   Optim.optimize(cost, X0c, Optim.BFGS())
  # else
  #   Optim.optimize(cost, X0c, Optim.NelderMead())
  # end
  alg = islen1 ? Optim.BFGS() : Optim.NelderMead()
  X0 = hat(M, ϵ, X0c)
  p0 = exp(M, ϵ, X0)
  r = Optim.optimize(Xc -> cost(p0, X0, Xc), X0c, alg)
  if !Optim.converged(r)
    @debug "Optim did not converge:" r
  end
  return exp(M, ϵ, hat(M, ϵ, r.minimizer))
end

#TODO Consolidate with _solveLambdaNumeric, see #1374
function _solveLambdaNumericMeas(
  fcttype::Union{F, <:Mixture{N_, F, S, T}},
  objResX::Function,
  residual::AbstractVector{<:Real},
  u0,#::AbstractVector{<:Real},
  variableType::InferenceVariable,
  islen1::Bool = false,
) where {N_, F <: AbstractManifoldMinimize, S, T}
  #
  # Assume measurement is on the tangent
  M = getManifold(variableType)#fcttype.M
  # the variable is a manifold point, we are working on the tangent plane in optim for now.
  ϵ = getPointIdentity(variableType)
  X0c = vee(M, ϵ, u0)

  function cost(X, Xc)
    hat!(M, X, ϵ, Xc)
    residual = objResX(X)
    return sum(residual .^ 2)
  end

  alg = islen1 ? Optim.BFGS() : Optim.NelderMead()
  X0 = hat(M, ϵ, X0c)
  r = Optim.optimize(Xc -> cost(X0, Xc), X0c, alg)
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
  # measurement_, # deprecate
  varParams,
  activehypo,
)
  #
  # FIXME, make thread safe (cache)
  # activevariables = view(ccwl.fullvariables, activehypo)
  activevariables = ccwl.fullvariables[activehypo]

  solveforidx = findfirst(==(ccwl.varidx), activehypo)

  return CalcFactor(
    _getusrfnc(ccwl),
    smpid,
    varParams,
    true,
    ccwl.dummyCache,
    tuple(activevariables...),
    solveforidx,
    getManifold(ccwl)
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
  target = view(ccwl.varValsAll[ccwl.varidx][smpid], ccwl.partialDims),
  measurement_ = ccwl.measurement;
  # fmd_::FactorMetadata = cpt_.factormetadata;
  _slack = nothing,
)
  #

  # TODO from obsolete _view:
  # Should be replaced with ccw.hypoParams::Tuple(hypo1, hypo2,...), made at construction and allows direct hypo lookup
  # DevNotes, also see new `hyporecipe` approach (towards consolidation CCW CPT FMd CF...)

  # build a view to the decision variable memory
  varValsHypo = ccwl.varValsAll[ccwl.activehypo]
  # tup = tuple(varParams...)
  # nms = keys(ccwl.varValsAll)[cpt_.activehypo]
  # varValsHypo = NamedTuple{nms,typeof(tup)}(tup)

  # prepare fmd according to hypo selection
  # FIXME must refactor (memory waste) and consolidate with CCW CPT FMd CF
  # FIXME move up out of smpid loop and only update bare minimal fields
  # _fmd_ = FactorMetadata( view(fmd_.fullvariables, cpt_.activehypo), 
  #                         view(fmd_.variablelist, cpt_.activehypo),
  #                         varValsHypo, #varParams, # view(fmd_.arrRef, cpt_.activehypo),
  #                         fmd_.solvefor,
  #                         fmd_.cachedata  )
  #
  # get the operational CalcFactor object
  cf = _buildCalcFactor(ccwl, smpid, varValsHypo, ccwl.activehypo)
  # new dev work on CalcFactor
  # cf = CalcFactor(ccwl.usrfnc!, _fmd_, smpid, 
  #                 varValsHypo)
  #

  # reset the residual vector
  fill!(ccwl.res, 0.0) # Roots->xDim | Minimize->zDim

  # build static lambda
  unrollHypo! = if _slack === nothing
    () -> cf(measurement_[smpid], map(vvh -> getindex(vvh, smpid), varValsHypo)...)
  else
    # slack is used to shift the residual away from the natural "zero" tension position of a factor, 
    # this is useful when calculating factor gradients at a variety of param locations resulting in "non-zero slack" of the residual.
    # see `IIF.calcFactorResidualTemporary`
    # NOTE this minus operation assumes _slack is either coordinate or tangent vector element (not a manifold or group element)
    () ->
      cf(measurement_[smpid], map(vvh -> getindex(vvh, smpid), varValsHypo)...) .- _slack
  end

  return unrollHypo!, target
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
  - perturb is necessary for NLsolve cases, and smaller than 1e-10 will result in test failure
- Also incorporates the active hypo lookup

DevNotes
- TODO testshuffle is now obsolete, should be removed
- TODO perhaps consolidate perturbation with inflation or nullhypo
"""
function _solveCCWNumeric!(
  ccwl::Union{CommonConvWrapper{F}, CommonConvWrapper{Mixture{N_, F, S, T}}};
  perturb::Real = 1e-10,
  # testshuffle::Bool=false,
  _slack = nothing,
) where {N_, F <: AbstractRelative, S, T}
  #

  #
  # thrid = Threads.threadid()

  smpid = ccwl.particleidx
  # cannot Nelder-Mead on 1dim, partial can be 1dim or more but being conservative.
  islen1 = length(ccwl.partialDims) == 1 || ccwl.partial
  # islen1 = length(cpt_.X[:, smpid]) == 1 || ccwl.partial

  # build the pre-objective function for this sample's hypothesis selection
  unrollHypo!, target = _buildCalcFactorLambdaSample(ccwl, smpid; _slack = _slack)

  # broadcast updates original view memory location
  ## using CalcFactor legacy path inside (::CalcFactor)
  _hypoObj = (x) -> (target .= x; unrollHypo!())

  # TODO small off-manifold perturbation is a numerical workaround only, make on-manifold requires RoME.jl #244
  # use all element dimensions : ==> 1:ccwl.xDim
  target .+= _perturbIfNecessary(getFactorType(ccwl), length(target), perturb)

  sfidx = ccwl.varidx
  # do the parameter search over defined decision variables using Minimization
  X = ccwl.varValsAll[sfidx][smpid][ccwl.partialDims]
  retval = _solveLambdaNumeric(
    getFactorType(ccwl), 
    _hypoObj, 
    ccwl.res, 
    X, 
    islen1
  )

  # Check for NaNs
  if sum(isnan.(retval)) != 0
    @error "$(ccwl.usrfnc!), ccw.thrid_=$(thrid), got NaN, smpid = $(smpid), r=$(retval)\n"
    return nothing
  end

  # insert result back at the correct variable element location
  ccwl.varValsAll[sfidx][smpid][ccwl.partialDims] .= retval

  return nothing
end
# brainstorming
# should only be calling a new arg list according to activehypo at start of particle
# Try calling an existing lambda
# sensitive to which hypo of course , see #1024
# need to shuffle content inside .cpt.fmd as well as .varValsAll accordingly
#

function _solveCCWNumeric!(
  ccwl::Union{CommonConvWrapper{F}, CommonConvWrapper{Mixture{N_, F, S, T}}};
  perturb::Real = 1e-10,
  # testshuffle::Bool=false,
  _slack = nothing,
) where {N_, F <: AbstractManifoldMinimize, S, T}
  #
  #   # FIXME, move this check higher and out of smpid loop
  # _checkErrorCCWNumerics(ccwl, testshuffle)

  #
  thrid = Threads.threadid()

  smpid = ccwl.particleidx
  # cannot Nelder-Mead on 1dim, partial can be 1dim or more but being conservative.
  islen1 = length(ccwl.partialDims) == 1 || ccwl.partial

  # build the pre-objective function for this sample's hypothesis selection
  unrollHypo!, target = _buildCalcFactorLambdaSample(
    ccwl,
    smpid,
    view(ccwl.varValsAll[ccwl.varidx], smpid);
    _slack = _slack,
  )

  # broadcast updates original view memory location
  ## using CalcFactor legacy path inside (::CalcFactor)
  # _hypoObj = (x) -> (target.=x; unrollHypo!())
  function _hypoObj(x)
    target[] .= x
    return unrollHypo!()
  end

  # TODO small off-manifold perturbation is a numerical workaround only, make on-manifold requires RoME.jl #244
  # use all element dimensions : ==> 1:ccwl.xDim
  # F <: AbstractRelativeRoots && (target .+= _perturbIfNecessary(getFactorType(ccwl), length(target), perturb))

  # do the parameter search over defined decision variables using Minimization
  sfidx = ccwl.varidx
  X = ccwl.varValsAll[sfidx][smpid]
  retval = _solveLambdaNumeric(
    getFactorType(ccwl),
    _hypoObj,
    ccwl.res,
    X,
    getVariableType(ccwl.fullvariables[sfidx]), # ccwl.vartypes[sfidx](), # only used for getting variable manifold and identity_element
    islen1,
  )

  # Check for NaNs
  # FIXME isnan for manifold ProductRepr
  # if sum(isnan.(retval)) != 0
  # @error "$(ccwl.usrfnc!), ccw.thrid_=$(thrid), got NaN, smpid = $(smpid), r=$(retval)\n"
  # return nothing
  # end

  # FIXME insert result back at the correct variable element location
  ccwl.varValsAll[sfidx][smpid] .= retval

  return nothing
end

#
