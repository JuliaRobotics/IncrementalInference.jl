## ================================================================================================
## Manifolds.jl Consolidation
## TODO: Still to be completed and tested.
## ================================================================================================
# struct ManifoldsVector <: Optim.Manifold
#   manis::Vector{Manifold}
# end

# Base.getindex(mv::ManifoldsVector, inds...) = getindex(mv.mani, inds...)
# Base.setindex!(mv, X, inds...) =  setindex!(mv.mani, X, inds...)

# function ManifoldsVector(fg::AbstractDFG, varIds::Vector{Symbol})
#   manis = Bool[]
#   for k = varIds
#     push!(manis, getVariableType(fg, k) |> getManifold)
#   end
#   ManifoldsVector(manis)
# end

# function Optim.retract!(manis::ManifoldsVector, x)
#   for (i,M) = enumerate(manis)
#     x[i] = project(M, x[i])
#   end
#   return x 
# end
# function Optim.project_tangent!(manis::ManifoldsVector, G, x)
#   for (i, M) = enumerate(manis)
#     G[i] = project(M, x[i], G)
#   end
#   return G
# end


##==============================================================================
## Old parametric kept for comparason until code is stabilized
##==============================================================================

"""
    $SIGNATURES

Batch solve a Gaussian factor graph using Optim.jl. Parameters can be passed directly to optim.
Notes:
  - Only :Euclid and :Circular manifolds are currently supported, own manifold are supported with `algorithmkwargs` (code may need updating though)
"""
function solveGraphParametric2(
  fg::AbstractDFG;
  computeCovariance::Bool = true,
  solvekey::Symbol = :parametric,
  autodiff = :forward,
  algorithm = Optim.BFGS,
  algorithmkwargs = (), # add manifold to overwrite computed one
  options = Optim.Options(;
    allow_f_increases = true,
    time_limit = 100,
    # show_trace = true,
    # show_every = 1,
  ),
)

  #Other options
  # options = Optim.Options(time_limit = 100,
  #                     iterations = 1000,
  #                     show_trace = true,
  #                     show_every = 1,
  #                     allow_f_increases=true,
  #                     g_tol = 1e-6,
  #                     )
  # Example for useing Optim's manifold functions
  # mc_mani = IIF.MixedCircular(fg, varIds)
  # alg = algorithm(;manifold=mc_mani, algorithmkwargs...)

  varIds = listVariables(fg)

  flatvar = FlatVariables(fg, varIds)

  for vId in varIds
    p = getVariableSolverData(fg, vId, solvekey).val[1]
    flatvar[vId] = getCoordinates(getVariableType(fg, vId), p)
  end

  initValues = flatvar.X
  # initValues .+= randn(length(initValues))*0.0001

  alg = algorithm(; algorithmkwargs...)

  cfd = calcFactorMahalanobisDict(fg)
  tdtotalCost = Optim.TwiceDifferentiable(
    (x) -> _totalCost(fg, cfd, flatvar, x),
    initValues;
    autodiff = autodiff,
  )

  result = Optim.optimize(tdtotalCost, initValues, alg, options)
  rv = Optim.minimizer(result)

  Σ = if computeCovariance
    H = Optim.hessian!(tdtotalCost, rv)
    pinv(H)
  else
    N = length(initValues)
    zeros(N, N)
  end

  d = Dict{Symbol, NamedTuple{(:val, :cov), Tuple{Vector{Float64}, Matrix{Float64}}}}()

  for key in varIds
    r = flatvar.idx[key]
    push!(d, key => (val = rv[r], cov = Σ[r, r]))
  end

  return d, result, flatvar.idx, Σ
end


##==============================================================================
## Deprecate code below before v0.34
##==============================================================================

# function CommonConvWrapper(
#   usrfnc::T,
#   fullvariables, #::Tuple ::Vector{<:DFGVariable};
#   varValsAll::Tuple,
#   X::AbstractVector{P}; #TODO remove X completely
#   # xDim::Int = size(X, 1),
#   userCache::CT = nothing,
#   manifold = getManifold(usrfnc),
#   partialDims::AbstractVector{<:Integer} = 1:length(X),
#   partial::Bool = false,
#   nullhypo::Real = 0,
#   inflation::Real = 3.0,
#   hypotheses::H = nothing,
#   certainhypo = nothing,
#   activehypo = collect(1:length(varValsAll)),
#   measurement::AbstractVector = Vector(Vector{Float64}()),
#   varidx::Int = 1,
#   particleidx::Int = 1,
#   res::AbstractVector{<:Real} = zeros(manifold_dimension(manifold)), # zDim
#   gradients = nothing,
# ) where {T <: AbstractFactor, P, H, CT}
#   #
#   return CommonConvWrapper(
#     usrfnc,
#     tuple(fullvariables...),
#     varValsAll,
#     userCache,
#     manifold,
#     partialDims,
#     partial,
#     # xDim,
#     float(nullhypo),
#     float(inflation),
#     hypotheses,
#     certainhypo,
#     activehypo,
#     measurement,
#     Ref(varidx),
#     Ref(particleidx),
#     res,
#     gradients,
#   )
# end

# function approxConvOnElements!(
#   ccwl::Union{CommonConvWrapper{F}, CommonConvWrapper{Mixture{N_, F, S, T}}},
#   elements::Union{Vector{Int}, UnitRange{Int}},
#   ::Type{<:MultiThreaded},
#   _slack = nothing,
# ) where {N_, F <: AbstractRelative, S, T}
#   #
#   return error(
#     "MultiThreaded `approxConvOnElements!` is deprecated and will soon be replaced",
#   )
#   # Threads.@threads for n in elements
#   #   # ccwl.thrid_ = Threads.threadid()
#   #   ccwl.cpt[Threads.threadid()].particleidx = n

#   #   # ccall(:jl_, Nothing, (Any,), "starting loop, thrid_=$(Threads.threadid()), partidx=$(ccwl.cpt[Threads.threadid()].particleidx)")
#   #   _solveCCWNumeric!( ccwl, _slack=_slack)
#   # end
#   # nothing
# end

# function approxConvOnElements!(
#   ccwl::Union{CommonConvWrapper{F}, CommonConvWrapper{Mixture{N_, F, S, T}}},
#   elements::Union{Vector{Int}, UnitRange{Int}},
#   _slack = nothing,
# ) where {N_, F <: AbstractRelative, S, T}
#   #
#   return approxConvOnElements!(ccwl, elements, ccwl.threadmodel, _slack)
# end

# more legacy, dont delete yet
function Base.getproperty(ccw::CommonConvWrapper, f::Symbol)
  if f == :threadmodel
    error("CommonConvWrapper.threadmodel is obsolete")
    # return SingleThreaded
  elseif f == :params
    error("CommonConvWrapper.params is deprecated, use .varValsAll instead")
    return ccw.varValsAll
  elseif f == :vartypes
    @warn "CommonConvWrapper.vartypes is deprecated, use typeof.(getVariableType.(ccw.fullvariables) instead" maxlog=3
    return typeof.(getVariableType.(ccw.fullvariables))
  else
    return getfield(ccw, f)
  end
end

##==============================================================================
## Deprecate code below before v0.35
##==============================================================================

##