
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

function _solveLambdaNumeric(
  fcttype::Union{F, <:Mixture{N_, F, S, T}},
  objResX::Function,
  residual::AbstractVector{<:Real},
  u0,#::AbstractVector{<:Real},
  variableType::InferenceVariable,
  islen1::Bool = false,
) where {N_, F <: AbstractManifoldMinimize, S, T}

  return _solveCCWNumeric_test_SA(fcttype, objResX, residual, u0, variableType, islen1)
  # return _solveLambdaNumeric_test_optim_manifold(fcttype, objResX, residual, u0, variableType, islen1)

end

function _solveLambdaNumeric_original(
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

# 1.355700 seconds (11.78 M allocations: 557.677 MiB, 6.96% gc time)
function _solveCCWNumeric_test_SA(
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

  X0c = zero(MVector{getDimension(M),Float64})
  X0c .= vee(M, u0, log(M, ϵ, u0))

  #TODO check performance
  function cost(Xc)
    X = hat(M, ϵ, Xc)
    p = exp(M, ϵ, X)  
    residual = objResX(p)
    # return sum(residual .^ 2)
    return sum(abs2, residual) #TODO maybe move this to CalcFactorNormSq
  end

  alg = islen1 ? Optim.BFGS() : Optim.NelderMead()

  r = Optim.optimize(cost, X0c, alg)
  if !Optim.converged(r)
    # TODO find good way for a solve to store diagnostics about number of failed converges etc.
    @warn "Optim did not converge (maxlog=10):" r maxlog=10
  end
  return exp(M, ϵ, hat(M, ϵ, r.minimizer))
end

# sloooowwww and does not always converge, unusable slow with gradient
# NelderMead 5.513693 seconds (38.60 M allocations: 1.613 GiB, 6.62% gc time)
function _solveLambdaNumeric_test_optim_manifold(
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
  ϵ = getPointIdentity(variableType)

  function cost(p)
    residual = objResX(p)
    return sum(residual .^ 2)
  end

  alg = islen1 ? Optim.BFGS(;manifold=ManifoldWrapper(M)) : Optim.NelderMead(;manifold=ManifoldWrapper(M))
  # alg = Optim.ConjugateGradient(; manifold=ManifoldWrapper(M))
  # alg = Optim.BFGS(; manifold=ManifoldWrapper(M))

  # r_backend = ManifoldDiff.TangentDiffBackend(
  #   ManifoldDiff.FiniteDifferencesBackend()
  # )
  
  # ## finitediff gradient
  # function costgrad_FD!(X,p)
  #   copyto!(X, ManifoldDiff.gradient(M, cost, p, r_backend))
  #   X
  # end

  u0_m = allocate(M, u0)
  u0_m .= u0 
  # r = Optim.optimize(cost, costgrad_FD!, u0_m, alg)
  r = Optim.optimize(cost, u0_m, alg)
  
  if !Optim.converged(r)
    @warn "Optim did not converge:" r maxlog=10
  end

  return r.minimizer
  # return exp(M, ϵ, hat(M, ϵ, r.minimizer))
end

#TODO Consolidate with _solveLambdaNumeric, see #1374
#TODO _solveLambdaNumericMeas assumes a measurement is always a tangent vector, confirm.
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

  function cost(Xc)
    X = hat(M, ϵ, Xc)
    residual = objResX(X)
    return sum(residual .^ 2)
  end

  alg = islen1 ? Optim.BFGS() : Optim.NelderMead()

  r = Optim.optimize(cost, X0c, alg)
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
)
  #
  # FIXME, make thread safe (cache)
  # activevariables = view(ccwl.fullvariables, activehypo)
  activevariables = ccwl.fullvariables[activehypo]

  solveforidx = findfirst(==(ccwl.varidx[]), activehypo)

  return CalcFactorNormSq(
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
  target, # partials no longer on coordinates at this level # = view(destVarVals[smpid], ccwl.partialDims), # target = view(ccwl.varValsAll[ccwl.varidx[]][smpid], ccwl.partialDims),
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

  _getindex_anyn(vec, n) = begin
    len = length(vec)
    # 1:len or any random element in that range
    getindex(vec, n <= len ? n : rand(1:len) )
  end

  # build static lambda
  unrollHypo! = if _slack === nothing
    # DESIGN DECISION WAS MADE THAT CALCFACTOR CALLS DO NOT DO INPLACE CHANGES TO ARGUMENTS, INSTEAD USING ISBITSTYPEs!!!!!!!!!
    () -> cf(measurement_[smpid], map(vvh -> _getindex_anyn(vvh, smpid), varValsHypo)...)
  else
    # slack is used to shift the residual away from the natural "zero" tension position of a factor, 
    # this is useful when calculating factor gradients at a variety of param locations resulting in "non-zero slack" of the residual.
    # see `IIF.calcFactorResidualTemporary`
    # NOTE this minus operation assumes _slack is either coordinate or tangent vector element (not a manifold or group element)
    () ->
      cf(measurement_[smpid], map(vvh -> _getindex_anyn(vvh, smpid), varValsHypo)...) .- _slack
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
  unrollHypo!, _ = _buildCalcFactorLambdaSample(
      # destVarVals,
      ccwl, 
      smpid, 
      target,
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
  # SUPER IMPORTANT ON PARTIALS, RESIDUAL FUNCTION MUST DEAL WITH PARTIAL AND WILL GET FULL VARIABLE POINTS REGARDLESS
  unrollHypo!, target = _buildCalcFactorLambdaSample(
    # destVarVals,
    ccwl,
    smpid,
    view(ccwl.varValsAll[][ccwl.varidx[]], smpid), # SUPER IMPORTANT, this `target` is mem pointer that will be updated by optim library
    ccwl.measurement,
    _slack,
  )

  # broadcast updates original view memory location
  ## using CalcFactor legacy path inside (::CalcFactor)
  # _hypoObj = (x) -> (target.=x; unrollHypo!())
  function _hypoObj(x)
    target[] = x
    return unrollHypo!()
  end

  # TODO small off-manifold perturbation is a numerical workaround only, make on-manifold requires RoME.jl #244
  # _perturbIfNecessary(getFactorType(ccwl), length(target), perturb)

  # do the parameter search over defined decision variables using Minimization
  sfidx = ccwl.varidx[]
  X = ccwl.varValsAll[][ccwl.varidx[]][smpid]
  retval = _solveLambdaNumeric(
    getFactorType(ccwl),
    _hypoObj,
    ccwl.res,
    X,
    getVariableType(ccwl.fullvariables[sfidx]), # ccwl.vartypes[sfidx](), # only used for getting variable manifold and identity_element
    islen1,
  )

  # TBD Check for NaNs

  # NOTE insert result back at the correct variable element location
  ccwl.varValsAll[][ccwl.varidx[]][smpid] = retval

  return nothing
end

#
