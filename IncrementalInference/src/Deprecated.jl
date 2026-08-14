
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

## ================================================================================================
## Deprecated in v0.37 -- the Optim.jl based parametric solve
##
## Superseded by the Manopt Riemannian Levenberg-Marquardt path:
##
##   batch      `solveGraphParametric!`   (`solve_RLM`, src/parametric/services/ParametricManopt.jl)
##   conditional `solve_RLM_conditional` / `solve_RLM_marginal` / `solve_RLM_propagate`
##   init        `autoinitParametric!`
##   tree        `solveTree!(fg; algorithm = :parametric)`
##
## ================================================================================================

# ================================================================================================
# FlatVariables - used for packing variables for optimization
# ================================================================================================

struct FlatVariables{T <: Real}
  X::Vector{T}
  idx::OrderedDict{Symbol, UnitRange{Int}}
end

function FlatVariables(fg::AbstractDFG, varIds::Vector{Symbol})
  index = 1
  idx = OrderedDict{Symbol, UnitRange{Int}}()
  for vid in varIds
    v = getVariable(fg, vid)
    dims = getDimension(v)
    idx[vid] = index:(index + dims - 1)
    index += dims
  end
  return FlatVariables(Vector{Float64}(undef, index - 1), idx)
end

function Base.setindex!(
  flatVar::FlatVariables{T},
  val::AbstractVector{T},
  vId::Symbol,
) where {T <: Real}
  if length(val) == length(flatVar.idx[vId])
    flatVar.X[flatVar.idx[vId]] .= val
  else
    error("array could not be broadcast to match destination")
  end
end

function Base.getindex(flatVar::FlatVariables{T}, vId::Symbol) where {T <: Real}
  return flatVar.X[flatVar.idx[vId]]
end

# ================================================================================================
# Parametric solve with Mahalanobis distance - CalcFactor
# ================================================================================================

function CalcFactorMahalanobis(fg, fct::FactorCompute)
  fac_func = getObservation(fct)
  varOrder = collect(getVariableOrder(fct))

  # NOTE, use getMeasurementParametric on FactorCompute{<:CCW} to allow special cases like OAS factors
  _meas, _iΣ = getFactorMeasurementParametric(fct) # fac_func
  
  # make sure its a tuple TODO Fix with mixture rework #1504
  meas = typeof(_meas) <: Tuple ? _meas : (_meas,)
  iΣ = typeof(_iΣ) <: Tuple ? _iΣ : (_iΣ,)

  cache = preambleCache(fg, getVariable.(fg, varOrder), getObservation(fct))

  multihypo = fct.hyper.multihypo
  nullhypo = fct.hyper.nullhypo

  # FIXME, type instability
  if length(multihypo) > 0
    special = MaxMultihypo(multihypo)
  elseif nullhypo > 0
    special = MaxNullhypo(nullhypo)
  elseif fac_func isa Mixture
    special = MaxMixture(fac_func.diversity.p, Ref(0))
  else
    special = nothing
  end

  return CalcFactorMahalanobis(fct.label, fac_func, cache, varOrder, meas, iΣ, special)
end

# This is where the actual parametric calculation happens, CalcFactor equivalent for parametric
# function (cfp::CalcFactorMahalanobis{FT, 1, C, MEAS, D, L, Nothing})(variables...) where {FT, C, MEAS, D, L, Nothing}# AbstractArray{T} where T <: Real
#   # call the user function
#   res = cfp.calcfactor!(cfp.meas..., variables...)
#   # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  # k = dim(μ)
#   return res' * cfp.iΣ[1] * res
# end

# function (cfm::CalcFactorMahalanobis)(variables...)
#   meas = cfm.meas
#   points = map(idx->p[idx], cfm.varOrderIdxs)
#   return cfm.sqrt_iΣ * cfm(meas, points...)
# end

function calcFactorMahalanobisDict(fg)
  calcFactors = OrderedDict{Symbol, CalcFactorMahalanobis}()
  for fct in getFactors(fg)
    # skip non-numeric prior
    getObservation(fct) isa MetaPrior ? continue : nothing
    calcFactors[fct.label] = CalcFactorMahalanobis(fg, fct)
  end
  return calcFactors
end

function calcFactorMahalanobisVec(fg)
  factypes, typedict, alltypes = getFactorTypesCount(getFactors(fg))
  
  # skip non-numeric prior (MetaPrior)
  #TODO test... remove MetaPrior{T} something like this
  metaPriorKeys = filter(k->contains(string(k), "MetaPrior"), collect(keys(alltypes)))
  delete!.(Ref(alltypes), metaPriorKeys)

  parts = map(values(alltypes)) do labels
    map(getFactor.(fg, labels)) do fct
      CalcFactorMahalanobis(fg, fct)
    end
  end
  parts_tuple = (parts...,)
  return ArrayPartition{CalcFactorMahalanobis, typeof(parts_tuple)}(parts_tuple)
end

struct GraphSolveBuffers{T <: Real, U}
  ϵ::U
  p::U
  X::U
  Xc::Vector{T}
end

function GraphSolveBuffers(@nospecialize(M), ::Type{T}) where {T}
  ϵ = getPointIdentity(M, T)
  p = deepcopy(ϵ)# allocate_result(M, getPointIdentity)
  X = deepcopy(ϵ) #allcoate(p)
  #FIXME update to ProductLieGroup first, but only 2 groups supported.
  Xc = Manifolds.get_coordinates(M, ϵ, X, DefaultOrthogonalBasis())
  # Xc = vee(LieGroup(M), X)
  return GraphSolveBuffers(ϵ, p, X, Xc)
end

struct GraphSolveContainer{CFT}
  M::AbstractManifold # ProductManifold or ProductGroup
  buffers::OrderedDict{DataType, GraphSolveBuffers}
  varTypes::Vector{DataType}
  varTypesIds::OrderedDict{DataType, Vector{Symbol}}
  varOrderDict::OrderedDict{Symbol, Tuple{Int, Vararg{Int}}}
  cfv::ArrayPartition{CalcFactorMahalanobis, CFT}
end

function GraphSolveContainer(fg)
  M, varTypes, varTypesIds = buildGraphSolveManifold(fg)
  varTypesIndexes = ArrayPartition(values(varTypesIds)...)
  buffs = OrderedDict{DataType, GraphSolveBuffers}()
  cfvec = calcFactorMahalanobisVec(fg)

  varOrderDict = OrderedDict{Symbol, Tuple{Int, Vararg{Int}}}()
  for cfp in cfvec
    fid = cfp.faclbl
    varOrder = cfp.varOrder
    var_idx = map(varOrder) do v
      return findfirst(==(v), varTypesIndexes)
    end
    varOrderDict[fid] = tuple(var_idx...)
  end

  return GraphSolveContainer(M, buffs, varTypes, varTypesIds,  varOrderDict, cfvec)
end

function getGraphSolveCache!(gsc::GraphSolveContainer, ::Type{T}) where {T <: Real}
  cache = gsc.buffers
  M = gsc.M
  val = get!(cache, T) do
    @debug "cache miss, cacheing" T
    return GraphSolveBuffers(M, T)
  end
  return val
end

function _toPoints2!(
  M::AbstractManifold,
  buffs::GraphSolveBuffers{T, U},
  Xc::Vector{T},
) where {T, U}
  ϵ = buffs.ϵ
  p = buffs.p
  X = buffs.X
  get_vector!(M, X, ϵ, Xc, DefaultOrthogonalBasis())
  exp!(M, p, ϵ, X)
  return p::U
end

function cost_cfp(
  cfp::CalcFactorMahalanobis,
  p::AbstractArray{T},
  vi::NTuple{N, Int},
) where {T,N}
  # cfp(map(v->p[v],vi)...)
  res = cfp(cfp.meas..., map(v->p[v],vi)...)
  # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  # k = dim(μ)
  return res' * cfp.iΣ[1] * res

end
# function cost_cfp(
#   @nospecialize(cfp::CalcFactorMahalanobis),
#   @nospecialize(p::AbstractArray),
#   vi::NTuple{1, Int},
# )
#   return cfp(p[vi[1]])
# end
# function cost_cfp(
#   @nospecialize(cfp::CalcFactorMahalanobis),
#   @nospecialize(p::AbstractArray),
#   vi::NTuple{2, Int},
# )
#   return cfp(p[vi[1]], p[vi[2]])
# end
# function cost_cfp(
#   @nospecialize(cfp::CalcFactorMahalanobis),
#   @nospecialize(p::AbstractArray),
#   vi::NTuple{3, Int},
# )
#   return cfp(p[vi[1]], p[vi[2]], p[vi[3]])
# end


# function (gsc::GraphSolveContainer)(f::Vector{T}, Xc::Vector{T}, ::Val{true}) where T <: Real
#   #
#   buffs = getGraphSolveCache!(gsc, T)

#   cfdict = gsc.cfdict
#   varOrderDict = gsc.varOrderDict

#   M = gsc.M 

#   p = _toPoints2!(M, buffs, Xc)

#   for (i,(fid, cfp)) in enumerate(cfdict)
#     varOrder_idx = varOrderDict[fid]

#     # call the user function
#     f[i] = cost_cfp(cfp, p, varOrder_idx)/2
#   end

#   return f
# end

# the cost function
function (gsc::GraphSolveContainer)(Xc::Vector{T}) where {T <: Real}
  #
  buffs = getGraphSolveCache!(gsc, T)

  varOrderDict = gsc.varOrderDict

  M = gsc.M

  p = _toPoints2!(M, buffs, Xc)
  
  obj = mapreduce(+, eachindex(gsc.cfv)) do i
    cfp = gsc.cfv[i]
    varOrder_idx = varOrderDict[cfp.faclbl]
    # # call the user function
    cost::T = cost_cfp(cfp, p, varOrder_idx)
    
    return cost
  end

  return obj / 2
end


# FIXME, deprecate and improve legacy use of `MultiThreaded` type
struct MultiThreaded end

function (gsc::GraphSolveContainer)(Xc::Vector{T}, ::MultiThreaded) where {T <: Real}
  #
  buffs = getGraphSolveCache!(gsc, T)

  cfdict = gsc.cfdict
  varOrderDict = gsc.varOrderDict

  M = gsc.M

  p = _toPoints2!(M, buffs, Xc)

  #NOTE multi threaded option
  obj = zeros(T, (Threads.nthreads()))
  Threads.@threads for fid in collect(keys(cfdict))
    cfp = cfdict[fid]

    #NOTE single thread option
    # obj::T = zero(T)
    # for (fid, cfp) in cfdict 

    varOrder_idx = varOrderDict[fid]

    # call the user function
    retval = cost_cfp(cfp, p, varOrder_idx)

    #NOTE multi threaded option
    obj[Threads.threadid()] += retval
    # NOTE single thread option
    # obj += retval
  end

  # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  # k = dim(μ)

  #NOTE multi threaded option
  return sum(obj) / 2
  # NOTE single thread option
  # return obj/2
end

#fg = generateCanonicalFG_Honeycomb!()

# copy variables from graph
function initPoints!(p, gsc, fg::AbstractDFG, solveKey = :parametric)
  for (i, vartype) in enumerate(gsc.varTypes)
    varIds = gsc.varTypesIds[vartype]
    for (j, vId) in enumerate(varIds)
      p[gsc.M, i][j] = DFG.refMeans(getState(fg, vId, solveKey))[1]
    end
  end
end

"""
    $SIGNATURES
Add parametric solver to fg, batch solve using [`solveGraphParametric`](@ref) and update fg.
"""
function solveGraphParametricOptim!(
  fg::AbstractDFG;
  init::Bool = true,
  solveKey::Symbol = :parametric,
  initSolveKey::Union{Nothing, Symbol} = nothing,
  verbose = false,
  kwargs...
)
  Base.depwarn(
    "`solveGraphParametricOptim!` (Optim.jl based) is deprecated, use `solveGraphParametric!` (Manopt RLM) or `solveTree!(fg; algorithm = :parametric)`.",
    :solveGraphParametricOptim!,
  )
  # make sure variables has solverData, see #1637
  makeSolverData!(fg; solveKey, parametric = true)
  if init
    if isnothing(initSolveKey)
      autoinitParametric!(fg; solveKey)
    else
      initParametricFrom!(fg, initSolveKey; parkey = solveKey)
    end
  end

  vardict, result, varIds, Σ = solveGraphParametricOptim(fg; solveKey, verbose, kwargs...)

  updateParametricSolution!(fg, vardict; solveKey)

  return vardict, result, varIds, Σ
end

function solveGraphParametricOptim(
  fg::AbstractDFG;
  verbose::Bool = false,
  computeCovariance::Bool = true,
  solveKey::Symbol = :parametric,
  autodiff = :forward,
  algorithm = Optim.BFGS,
  algorithmkwargs = (), # add manifold to overwrite computed one
  # algorithmkwargs = (linesearch=Optim.BackTracking(),), # add manifold to overwrite computed one
  options = Optim.Options(;
    allow_f_increases = true,
    time_limit = 100,
    # show_trace = true,
    # show_every = 1,
  ),
)
  Base.depwarn(
    "`solveGraphParametricOptim` (Optim.jl based) is deprecated, use `solveGraphParametric!` (Manopt RLM) or `solveTree!(fg; algorithm = :parametric)`.",
    :solveGraphParametricOptim,
  )
  # 
  # Build the container  
  gsc = GraphSolveContainer(fg)
  buffs = getGraphSolveCache!(gsc, Float64)

  M = gsc.M
  ϵ = buffs.ϵ
  p = buffs.p
  X = buffs.X
  Xc = buffs.Xc

  #initialize points in buffer from fg, TODO maybe do in constructor
  initPoints!(p, gsc, fg, solveKey)

  # log!(M, X, Identity(ProductOperation), p)
  # calculate initial coordinates vector for Optim
  log!(M, X, ϵ, p)
  #FIXME update to ProductLieGroup first, but only 2 groups supported.
  get_coordinates!(M, Xc, ϵ, X, DefaultOrthogonalBasis())
  # vee!(LieGroup(M), Xc, X)

  initValues = Xc
  #FIXME, for some reason we get NANs and adding a small random value works
  initValues .+= randn(length(Xc)) * 0.0001

  #optim setup and solve
  alg = algorithm(; algorithmkwargs...)

  tdtotalCost = Optim.TwiceDifferentiable(gsc, initValues; autodiff = autodiff)

  result = Optim.optimize(tdtotalCost, initValues, alg, options)
  !verbose ? nothing : @show(result)

  rv = Optim.minimizer(result)

  # optionally compute hessian for covariance
  Σ = if computeCovariance
    H = Optim.hessian!(tdtotalCost, rv)
    pinv(H)
  else
    N = length(initValues)
    zeros(N, N)
  end

  #TODO better return 

  #get point (p) values form results
  get_vector!(M, X, ϵ, rv, DefaultOrthogonalBasis())
  exp!(M, p, ϵ, X)

  #extract covariances from result
  # sigmas = getComponentsCovar(M, Σ)

  # d = OrderedDict{Symbol,NamedTuple{(:val, :cov),Tuple{Vector{Float64},Matrix{Float64}}}}()
  d = OrderedDict{Symbol, NamedTuple{(:val, :cov), Tuple{AbstractArray, Matrix{Float64}}}}()

  varIds = vcat(values(gsc.varTypesIds)...)
  varIdDict = FlatVariables(fg, varIds).idx
  for (i, key) in enumerate(varIds)
    r = varIdDict[key]
    push!(d, key => (val = p[i], cov = Σ[r, r]))
    # push!(d,key=>(val=p[i], cov=sigmas[i]))
  end

  return (opti = d, stat = result, varIds = varIdDict, Σ = Σ)
end
# Original
# ==============================

function _totalCost(fg, cfdict::OrderedDict{Symbol, <:CalcFactorMahalanobis}, flatvar, Xc)
  #
  obj = zero(eltype(Xc))
  for (fid, cfp) in cfdict
    varOrder = cfp.varOrder

    Xparams = [
      getPoint(getStateKind(fg, varId), view(Xc, flatvar.idx[varId])) for
      varId in varOrder
    ]

    # call the user function
    # retval = cfp(Xparams...)
    res = cfp(cfp.meas..., Xparams...)
    # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  # k = dim(μ)
    obj += 1 / 2 * res' * cfp.iΣ[1] * res
  end

  return obj
end

"""
$SIGNATURES
Solve for frontal values only with values in seprarators fixed
  
DevNotes
- WIP
- Relates to: https://github.com/JuliaRobotics/IncrementalInference.jl/issues/466#issuecomment-562556953
- Consolidation
  - Definitely with [`solveFactorParametric`](@ref)
  - Maybe with [`solveGraphParametric`](@ref)
    - https://github.com/JuliaRobotics/IncrementalInference.jl/pull/1588#issuecomment-1210406683
"""
function solveConditionalsParametric(
  fg::AbstractDFG,
  frontals::Vector{Symbol},
  separators::Vector{Symbol} = setdiff(listVariables(fg), frontals);
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
  Base.depwarn(
    "`solveConditionalsParametric` (Optim.jl based) is deprecated, use `solve_RLM_conditional` / `solve_RLM_marginal` / `solve_RLM_propagate`.",
    :solveConditionalsParametric,
  )
  varIds = [frontals; separators]

  sfg = issetequal(varIds, listVariables(fg)) ? fg : getSubgraph(fg, varIds, 1)

  flatvar = FlatVariables(fg, varIds)

  for vId in varIds
    p = DFG.refMeans(getState(fg, vId, solvekey))[1]
    flatvar[vId] = getCoordinates(getStateKind(fg, vId), p)
  end
  initValues = flatvar.X

  frontalsLength = sum(map(v -> getDimension(getVariable(fg, v)), frontals))

  # build variables for frontals and seperators
  # fX = view(initValues, 1:frontalsLength)
  fX = initValues[1:frontalsLength]
  # sX = view(initValues, (frontalsLength+1):length(initValues))
  sX = initValues[(frontalsLength + 1):end]

  alg = algorithm(; algorithmkwargs...)
  # alg = algorithm(; algorithmkwargs...)
  cfd = calcFactorMahalanobisDict(sfg)
  tdtotalCost = Optim.TwiceDifferentiable(
    (x) -> _totalCost(fg, cfd, flatvar, [x; sX]),
    fX;
    autodiff = autodiff,
  )

  # result = Optim.optimize((x)->_totalCost(fg, flatvar, [x;sX]), fX, alg, options)
  result = Optim.optimize(tdtotalCost, fX, alg, options)

  if !Optim.converged(result)
    @warn "Optim did not converge:" result maxlog=10
  end

  rv = Optim.minimizer(result)

  H = Optim.hessian!(tdtotalCost, rv)

  Σ = pinv(H)

  d = OrderedDict{Symbol, NamedTuple{(:val, :cov), Tuple{AbstractArray, Matrix{Float64}}}}()

  for key in frontals
    r = flatvar.idx[key]
    p = getPoint(getStateKind(fg, key), rv[r])
    push!(d, key => (val = p, cov = Σ[r, r]))
  end

  return (opti = d, stat = result, varIds = flatvar.idx, Σ = Σ)
end
# ================================================================================================
# UNDER DEVELOPMENT Parametric solveTree utils
# ================================================================================================

"""
    $SIGNATURES
Get the indexes for labels in FlatVariables
"""
function collectIdx(varinds, labels)
  idx = Int[]
  for lbl in labels
    append!(idx, varinds[lbl])
  end
  return idx
end

"""
    $SIGNATURES
Calculate the marginal distribution for a clique over subsetVarIds.
#FIXME update to support manifolds
"""
function calculateMarginalCliqueLikelihood(vardict, Σ, varindxs, subsetVarIds)
  μₘ = Float64[]
  for lbl in subsetVarIds
    append!(μₘ, vardict[lbl].val)
  end

  Aidx = collectIdx(varindxs, subsetVarIds)
  Σₘ = Σ[Aidx, Aidx]

  return createMvNormal(μₘ, Σₘ)
end

function autoinitParametricOptim!(
  fg,
  clique_order = getInitOrderParametric(fg);
  reinit = false,
  algorithm = Optim.NelderMead,
  algorithmkwargs = (initial_simplex = Optim.AffineSimplexer(0.025, 0.1),),
  kwargs...
)
  Base.depwarn(
    "`autoinitParametricOptim!` (Optim.jl based) is deprecated, use `autoinitParametric!`.",
    :autoinitParametricOptim!,
  )
  @showprogress for cliq in clique_order
    for vIdx in cliq.frontals
      autoinitParametricOptim!(fg, vIdx; reinit, algorithm, algorithmkwargs, kwargs...)
    end
  end
  return nothing
end

function autoinitParametricOptim!(dfg::AbstractDFG, initme::Symbol; kwargs...)
  return autoinitParametricOptim!(dfg, getVariable(dfg, initme); kwargs...)
end

function autoinitParametricOptim!(
  dfg::AbstractDFG,
  xi::VariableCompute;
  solveKey = :parametric,
  reinit::Bool = false,
  kwargs...,
)
  #

  initme = getLabel(xi)
  vnd = getState(xi, solveKey)
  # don't initialize a variable more than once
  if reinit || !isInitialized(xi, solveKey)

    # frontals - initme
    # separators - inifrom

    initfrom = ls2(dfg, initme)
    filter!(initfrom) do vl
      return isInitialized(dfg, vl, solveKey)
    end

    vardict, result, flatvars, Σ =
      solveConditionalsParametric(dfg, [initme], initfrom; kwargs...)

    val, cov = vardict[initme]

    updateSolverDataParametric!(vnd, val, cov)

    vnd.initialized = true

    # updateVariableSolverData!(dfg, xi, solveKey, true; warn_if_absent=false)    
    # updateVariableSolverData!(dfg, xi.label, getState(xi, solveKey), :graphinit, true, Symbol[]; warn_if_absent=false)
  else
    result = nothing
  end

  return result#isInitialized(xi, solveKey)
end

# ================================================================================================
# LazyCase based on LazyBufferCache from PreallocationTools.jl
# ================================================================================================

"""
  $SIGNATURES
A lazily allocated cache object.
"""
struct LazyCache{F <: Function}
  dict::Dict{Tuple{DataType, Symbol}, Any}
  fnc::F
end
function LazyCache(f::F = allocate) where {F <: Function}
  return LazyCache(Dict{Tuple{DataType, Symbol}, Any}(), f)
end

# override the [] method
function Base.getindex(cache::LazyCache, u::T, varname::Symbol) where {T}
  val = get!(cache.dict, (T, varname)) do
    return cache.fnc(u)
  end::T
  return val
end

function getCoordCache!(cache::LazyCache, M, T::DataType, varname::Symbol)
  val = get!(cache.dict, (T, varname)) do
    return Vector{T}(undef, manifold_dimension(M))
  end::Vector{T}
  return val
end