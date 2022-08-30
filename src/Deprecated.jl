
##==============================================================================
## LEGACY, towards Sidecar
##==============================================================================

"""
Converter: Prior -> PackedPrior::Dict{String, Any}

FIXME see DFG #590 for consolidation with Serialization and Marshaling
"""
function convert(::Type{Dict{String, Any}}, prior::IncrementalInference.Prior)
  @error("Obsolete, use pack/unpack converters instead")
  z = convert(Type{Dict{String, Any}}, prior.Z)
  return Packed_Factor([z], "Prior")
end

"""
Converter: PackedPrior::Dict{String, Any} -> Prior

FIXME see DFG #590 for consolidation on Serialization and Marshaling
"""
function convert(::Type{<:Prior}, prior::Dict{String, Any})
  @error("Obsolete, use pack/unpack converters instead")
  # Genericize to any packed type next.
  z = prior["measurement"][1]
  z = convert(DFG.getTypeFromSerializationModule(z["distType"]), z)
  return Prior(z)
end

##==============================================================================
## Deprecate code below before v0.32
##==============================================================================

# """
#     $SIGNATURES
# Get `.factormetadata` for each CPT in CCW for a specific factor in `fg`. 
# """
# _getFMdThread(ccw::CommonConvWrapper, 
#               thrid::Int=Threads.threadid()) = ccw.cpt[thrid].factormetadata
# #
# _getFMdThread(fc::Union{GenericFunctionNodeData,DFGFactor}, 
#               thrid::Int=Threads.threadid()) = _getFMdThread(_getCCW(fc), thrid)
# #
# _getFMdThread(dfg::AbstractDFG,
#               lbl::Symbol,
#               thrid::Int=Threads.threadid()) = _getFMdThread(_getCCW(dfg, lbl), thrid)
# #

##==============================================================================
## Deprecate code below before v0.31
##==============================================================================

@deprecate initManual!(w...; kw...) initVariable!(w...; kw...)

##==============================================================================
## Old parametric kept for comparason until code stabilize
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
