
# moved here from DistributedFactorGraphs.jl, replace with new way.
function typeModuleName(variableType::StateType)
    Base.depwarn("typeModuleName is obsolete", :typeModuleName)
    io = IOBuffer()
    ioc = IOContext(io, :module => DistributedFactorGraphs)
    show(ioc, typeof(variableType))
    return String(take!(io))
end

"""
    $(SIGNATURES)
Get a type from the serialization module.
"""
function getTypeFromSerializationModule(_typeString::AbstractString)
    @debug "DFG converting type string to Julia type" _typeString
    try
        # split the type at last `.`
        split_st = split(_typeString, r"\.(?!.*\.)")
        #if module is specified look for the module in main, otherwise use Main        
        if length(split_st) == 2
            m = getfield(Main, Symbol(split_st[1]))
        else
            m = Main
        end
        noparams = split(split_st[end], r"{")
        ret = if 1 < length(noparams)
            # fix #671, but does not work with specific module yet
            bidx = findfirst(r"{", split_st[end])[1]
            @error("getTypeFromSerializationModule eval obsolete")
            Core.eval(m, Base.Meta.parse("$(noparams[1])$(split_st[end][bidx:end])"))
        else
            getfield(m, Symbol(split_st[end]))
        end

        return ret

    catch ex
        @error "Unable to deserialize type $(_typeString)"
        io = IOBuffer()
        showerror(io, ex, catch_backtrace())
        err = String(take!(io))
        @error(err)
    end
    return nothing
end


## ================================================================================================
## Deprecated in v0.36
## ================================================================================================

#TODO this looks like dead code, should be removed
# TODO deprecate testshuffle
function _checkErrorCCWNumerics(
  ccwl::CommonConvWrapper{F},
  testshuffle::Bool = false,
) where {F <: AbstractRelativeMinimize}
  return nothing
end
function _checkErrorCCWNumerics(
  ccwl::CommonConvWrapper{F},
  testshuffle::Bool = false,
) where {F <: AbstractManifoldMinimize}
  return nothing
end

#TODO this looks like dead code, should be removed
function _perturbIfNecessary(
  fcttype::AbstractRelativeMinimize,
  len::Int = 1,
  perturbation::Real = 1e-10,
)
  return 0
end

function _perturbIfNecessary(
  fcttype::AbstractManifoldMinimize,
  len::Int = 1,
  perturbation::Real = 1e-10,
)
  return 0
end
#

# lets create all the vertices first and then deal with the elimination variables thereafter
function addBayesNetVerts!(dfg::AbstractDFG, elimOrder::Array{Symbol, 1})
  #
  for pId in elimOrder
    vert = DFG.getVariable(dfg, pId)
    if  getState(vert, :default).BayesNetVertID === nothing ||
        getState(vert, :default).BayesNetVertID == :_null # Special serialization case of nothing
      @debug "[AddBayesNetVerts] Assigning $pId.data.BayesNetVertID = $pId"
      getState(vert, :default).BayesNetVertID = pId
    else
      @warn "addBayesNetVerts -- Something is wrong, variable '$pId' should not have an existing Bayes net reference to '$(getState(vert, :default).BayesNetVertID)'"
    end
  end
end

## packing converters-----------------------------------------------------------
# heavy use of multiple dispatch for converting between packed and original data types during DB usage

# function convert(
#   ::Type{DFG.PackedFunctionNodeData{P}},
#   d::DFG.FunctionNodeData{T},
# ) where {P <: AbstractPackedObservation, T <: FactorCache}
#   error("TODO remove. PackedFunctionNodeData is obsolete")
#   return DFG.PackedFunctionNodeData(
#     d.eliminated,
#     d.potentialused,
#     d.edgeIDs,
#     convert(P, _getCCW(d).usrfnc!),
#     d.multihypo,
#     _getCCW(d).hyporecipe.certainhypo,
#     d.nullhypo,
#     d.solveInProgress,
#     d.inflation,
#   )  # extract two values from ccw for storage -- ccw thrown away
# end

## unpack converters------------------------------------------------------------
# see #1424
#TODO Consolidate: this looks alot like `getDefaultFactorData`
# function DFG.reconstFactorData(
#   dfg::AbstractDFG,
#   varOrder::AbstractVector{Symbol},
#   ::Type{<:DFG.GenericFunctionNodeData{<:CommonConvWrapper{F}}},
#   packed::DFG.GenericFunctionNodeData{<:AbstractPackedObservation},
# ) where {F <: AbstractObservation}

#   error("TODO remove. Obsolete: use `DFG.rebuildFactorCache!` and getDefaultFactorData instead.")
#   #
#   # TODO store threadmodel=MutliThreaded,SingleThreaded in persistence layer
#   usrfnc = convert(F, packed.fnc)
#   multihypo, nullhypo = parseusermultihypo(packed.multihypo, packed.nullhypo)

#   # IIF #1424
#   vars = map(f -> getVariable(dfg, f), varOrder)
#   userCache = preambleCache(dfg, vars, usrfnc)

#   # TODO -- improve _createCCW for hypotheses and certainhypo field recovery when deserializing
#   # reconstitute from stored data
#   # FIXME, add threadmodel=threadmodel
#   # FIXME https://github.com/JuliaRobotics/DistributedFactorGraphs.jl/issues/590#issuecomment-776838053
#   # FIXME dont know what manifolds to use in ccw
#   ccw = _createCCW(
#     vars,
#     usrfnc;
#     multihypo,
#     nullhypo,
#     certainhypo = packed.certainhypo,
#     inflation = packed.inflation,
#     userCache,
#     attemptGradients = getSolverParams(dfg).attemptGradients,
#     # Block recursion if NoSolverParams or if set to not attempt gradients.
#     _blockRecursion=
#       getSolverParams(dfg) isa NoSolverParams || 
#       !getSolverParams(dfg).attemptGradients,
#   )
#   #

#   # CommonConvWrapper{typeof(usrfnc)}
#   ret = DFG.FunctionNodeData{typeof(ccw)}(
#     packed.eliminated,
#     packed.potentialused,
#     packed.edgeIDs,
#     ccw,
#     packed.multihypo,
#     packed.certainhypo,
#     packed.nullhypo,
#     packed.solveInProgress,
#     packed.inflation,
#   )
#   #
#   return ret
# end

# function _getDimensionsPartial(data::DFG.GenericFunctionNodeData)
#   Base.depwarn(
#     "_getDimensionsPartial(data::GenericFunctionNodeData) is deprecated, use solvercache <: FactorCache instead",
#     :_getDimensionsPartial,
#   ) 
#   return _getCCW(data) |> _getDimensionsPartial
# end

# """
#     $SIGNATURES
# Get the CommonConvWrapper for this factor.
# """
# function _getCCW(gfnd::DFG.GenericFunctionNodeData)
#   error("_getCCW(gfnd::DFG.GenericFunctionNodeData) is deprecated, use DFG.getCache instead.")
# end

# _getZDim(fcd::DFG.GenericFunctionNodeData) = _getCCW(fcd) |> _getZDim
# DFG.getDimension(fct::DFG.GenericFunctionNodeData) = _getZDim(fct)

function sampleTangent(x::ManifoldKernelDensity, p = mean(x))
  error("sampleTangent(x::ManifoldKernelDensity, p) should be replaced by sampleTangent(M<:AbstractManifold, x::ManifoldKernelDensity, p)")
end

export setPPE!, setVariablePosteriorEstimates!
setPPE!(args...; kw...) = error("PPEs are obsolete (use `calcMeanMaxSuggested` provisionally), see DFG #1133")
setVariablePosteriorEstimates!(args...; kw...) = error("PPEs are obsolete (use `calcMeanMaxSuggested` provisionally), see DFG #1133")

@deprecate calcPPE(
  var::VariableCompute,
  varType::StateType = getStateKind(var);
  solveKey::Symbol = :default,
  kwargs...,
) calcMeanMaxSuggested(var, solveKey)

@deprecate calcPPE(
  dfg::AbstractDFG,
  label::Symbol;
  solveKey::Symbol = :default,
  kwargs...,
) calcMeanMaxSuggested(dfg, label, solveKey)

export calcVariablePPE
const calcVariablePPE = calcPPE

#FIXME The next functions use PPEs and should be updated or deprecated
# getPPESuggestedAll no external use
# findVariablesNear used in 1 rome example
"""
    $SIGNATURES

Return `::Tuple` with matching variable ID symbols and `Suggested` PPE values.

Related

getVariablePPE
"""
function getPPESuggestedAll(dfg::AbstractDFG, regexFilter::Union{Nothing, Regex} = nothing)
  #
  # get values
  vsyms = listVariables(dfg, regexFilter) |> sortDFG
  slamPPE = map(x -> getVariablePPE(dfg, x).suggested, vsyms)
  # sizes to convert to matrix
  rumax = zeros(Int, 2)
  for ppe in slamPPE
    rumax[2] = length(ppe)
    rumax[1] = maximum(rumax)
  end

  # populate with values
  XYT = zeros(length(slamPPE), rumax[1])
  for i = 1:length(slamPPE)
    XYT[i, 1:length(slamPPE[i])] = slamPPE[i]
  end
  return (vsyms, XYT)
end

"""
    $SIGNATURES

Find and return a `::Tuple` of variables and distances to `loc::Vector{<:Real}`.

Related

findVariablesNearTimestamp
"""
function findVariablesNear(
  dfg::AbstractDFG,
  loc::Vector{<:Real},
  regexFilter::Union{Nothing, Regex} = nothing;
  number::Int = 3,
)
  #

  xy = getPPESuggestedAll(dfg, regexFilter)
  dist = sum((xy[2][:, 1:length(loc)] .- loc') .^ 2; dims = 2) |> vec
  prm = (dist |> sortperm)[1:number]
  return (xy[1][prm], sqrt.(dist[prm]))
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
#     push!(manis, getStateKind(fg, k) |> getManifold)
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
    p = getState(fg, vId, solvekey).val[1]
    flatvar[vId] = getCoordinates(getStateKind(fg, vId), p)
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
#   # TODO get xDim = getDimension(getStateKind(Xi[sfidx])) but without having Xi
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
    @warn "CommonConvWrapper.vartypes is deprecated, use typeof.(getStateKind.(ccw.fullvariables) instead" maxlog=3
    return typeof.(getStateKind.(ccw.fullvariables))
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