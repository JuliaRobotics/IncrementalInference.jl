## ================================================================================================
## ================================================================================================

# TODO maybe upstream to DFG
DFG.MeanMaxPPE(solveKey::Symbol, suggested::SVector, max::SVector, mean::SVector) =
  DFG.MeanMaxPPE(solveKey, collect(suggested), collect(max), collect(mean))


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
## Deprecate code below before v0.37
##==============================================================================

@deprecate solveFactorParameteric(w...;kw...) solveFactorParametric(w...;kw...)

##==============================================================================
## Deprecate code below before v0.36
##==============================================================================

# function Base.isapprox(a::ProductRepr, b::ProductRepr; atol::Real = 1e-6)
#   #
#   for (i, a_) in enumerate(a.parts)
#     isapprox(a_, b.parts[i]; atol = atol) || (return false)
#   end
#   return true
# end

# exportimg(pl) = error("Please do `using Gadfly` to allow image export.")

# function _perturbIfNecessary(
#   fcttype::Union{F, <:Mixture{N_, F, S, T}},
#   len::Int = 1,
#   perturbation::Real = 1e-10,
# ) where {N_, F <: AbstractRelativeRoots, S, T}
#   return perturbation * randn(len)
# end

# function _checkErrorCCWNumerics(
#   ccwl::Union{CommonConvWrapper{F}, CommonConvWrapper{Mixture{N_, F, S, T}}},
#   testshuffle::Bool = false,
# ) where {N_, F <: AbstractRelativeRoots, S, T}
#   #
#   # error("<:AbstractRelativeRoots is obsolete, use one of the other <:AbstractRelative types instead.")
#   # TODO get xDim = getDimension(getVariableType(Xi[sfidx])) but without having Xi
#   if testshuffle || ccwl.partial
#     error(
#       "<:AbstractRelativeRoots factors with less or more measurement dimensions than variable dimensions have been discontinued, rather use <:AbstractManifoldMinimize.",
#     )
#   # elseif !(_getZDim(ccwl) >= ccwl.xDim && !ccwl.partial)
#   #   error("Unresolved numeric <:AbstractRelativeRoots solve case")
#   end
#   return nothing
# end

# function _solveLambdaNumeric(
#   fcttype::Union{F, <:Mixture{N_, F, S, T}},
#   objResX::Function,
#   residual::AbstractVector{<:Real},
#   u0::AbstractVector{<:Real},
#   islen1::Bool = false,
# ) where {N_, F <: AbstractRelativeRoots, S, T}
#   #

#   #
#   r = NLsolve.nlsolve((res, x) -> res .= objResX(x), u0; inplace = true) #, ftol=1e-14)

#   #
#   return r.zero
# end

# should probably deprecate the abstract type approach?
abstract type _AbstractThreadModel end

"""
$(TYPEDEF)
"""
struct SingleThreaded <: _AbstractThreadModel end
# """
# $(TYPEDEF)
# """
# struct MultiThreaded <: _AbstractThreadModel end


##==============================================================================
## Deprecate code below before v0.35
##==============================================================================


@deprecate _prepCCW(w...;kw...) _createCCW(w...;kw...)

predictbelief(w...;asPartial::Bool=false,kw...) = begin 
  @warn("predictbelief is deprecated, use propagateBelief instead")
  bel,ipc = propagateBelief(w...;asPartial,kw...)
  getPoints(bel), ipc
end


# more legacy, dont delete yet
function Base.getproperty(ccw::CommonConvWrapper, f::Symbol)
  if f == :threadmodel
    error("CommonConvWrapper.threadmodel is obsolete")
    # return SingleThreaded
  elseif f == :params
    error("CommonConvWrapper.params is deprecated, use .varValsAll instead")
    return ccw.varValsAll[]
  elseif f == :vartypes
    @warn "CommonConvWrapper.vartypes is deprecated, use typeof.(getVariableType.(ccw.fullvariables) instead" maxlog=3
    return typeof.(getVariableType.(ccw.fullvariables))
  elseif f == :hypotheses
    @warn "CommonConvWrapper.hypotheses is now under ccw.hyporecipe.hypotheses" maxlog=5
    return ccw.hyporecipe.hypotheses
  elseif f == :certainhypo
    @warn "CommonConvWrapper.certainhypo is now under ccw.hyporecipe.certainhypo" maxlog=5
    return ccw.hyporecipe.certainhypo
  elseif f == :activehypo
    @warn "CommonConvWrapper.activehypo is now under ccw.hyporecipe.activehypo" maxlog=5
    return ccw.hyporecipe.activehypo
  else
    return getfield(ccw, f)
  end
end


# function __init__()
#   # @require InteractiveUtils = "b77e0a4c-d291-57a0-90e8-8db25a27a240" include(
#   #   "services/RequireInteractiveUtils.jl",
#   # )
#   # @require Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004" include(
#   #   "services/EmbeddedPlottingUtils.jl",
#   # )
#   # @require DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa" include(
#   #   "ODE/DERelative.jl",
#   # )
#   # @require Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59" include(
#   #   "services/HeatmapSampler.jl",
#   # )

#   # # combining neural networks natively into the non-Gaussian  factor graph object
#   # @require Flux = "587475ba-b771-5e3f-ad9e-33799f191a9c" begin
#   #   # include("Flux/FluxModelsDistribution.jl")
#   #   include("Serialization/services/FluxModelsSerialization.jl") # uses BSON
#   # end
# end


##