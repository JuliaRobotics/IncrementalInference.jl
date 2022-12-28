
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

# more legacy, dont delete yet
function Base.getproperty(ccw::CommonConvWrapper, f::Symbol)
  if f == :threadmodel
    @warn "CommonConvWrapper.threadmodel is obsolete" maxlog=3
    return SingleThreaded
  elseif f == :params
    @warn "CommonConvWrapper.params is deprecated, use .varValsAll instead" maxlog=3
    return ccw.varValsAll
  elseif f == :vartypes
    @warn "CommonConvWrapper.vartypes is deprecated, use typeof.(getVariableType.(ccw.fullvariables) instead" maxlog=3
    return typeof.(getVariableType.(ccw.fullvariables))
  else
    return getfield(ccw, f)
  end
end


##==============================================================================
## Deprecate code below before v0.33
##==============================================================================

# export setThreadModel!
  # introduced for approximate convolution operations
export SingleThreaded, MultiThreaded

function setThreadModel!(fgl::AbstractDFG; model = IIF.SingleThreaded)
  #
  @error("Obsolete, ThreadModel types are no longer in use.")
  # for (key, id) in fgl.fIDs
  #   _getCCW(fgl, key).threadmodel = model
  # end
  return nothing
end

# should have been deleted in v0.31 but no harm in keeping this one a bit longer
@deprecate initManual!(w...; kw...) initVariable!(w...; kw...)

##==============================================================================
## Deprecate code below before v0.32
##==============================================================================

# NOTE Temporary fix? Overwrite deepcopy on MetaBayesTree to strip out copying the channels. 
# see https://github.com/JuliaRobotics/IncrementalInference.jl/issues/1530
# possible fix in https://github.com/JuliaLang/julia/pull/46406
import Base.deepcopy_internal
VERSION < v"1.8.1" && function Base.deepcopy_internal(bt::MetaBayesTree, stackdict::IdDict)
  if haskey(stackdict, bt)
    return stackdict[bt]
  end

  mg = bt.bt

  graph = deepcopy_internal(mg.graph, stackdict)
  vprops = deepcopy_internal(mg.vprops, stackdict)
  T = eltype(mg)
  # dropping all edge data
  eprops = Dict{MetaGraphs.SimpleEdge{T}, MetaGraphs.PropDict}()
  gprops = deepcopy_internal(mg.gprops, stackdict)
  weightfield = deepcopy_internal(mg.weightfield, stackdict)
  defaultweight = deepcopy_internal(mg.defaultweight, stackdict)
  metaindex = deepcopy_internal(mg.metaindex, stackdict)
  indices = deepcopy_internal(mg.indices, stackdict)

  mg_cpy = MetaDiGraph(
    graph,
    vprops,
    eprops,
    gprops,
    weightfield,
    defaultweight,
    metaindex,
    indices,
  )

  bt_cpy = MetaBayesTree(
    mg_cpy,
    bt.btid,
    deepcopy_internal(bt.frontals, stackdict),
    deepcopy_internal(bt.eliminationOrder, stackdict),
    bt.buildTime,
  )

  stackdict[bt] = bt_cpy
  return bt_cpy
end


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
