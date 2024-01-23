# ================================================================================================
## FlatVariables - used for packing variables for optimization
## ================================================================================================

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

## ================================================================================================
## Parametric Factors
## ================================================================================================

"""
    $SIGNATURES

Returns the parametric measurement for a factor as a tuple (measurement, inverse covariance) for parametric inference (assuming Gaussian).
Defaults to find the parametric measurement at field `Z`.
Notes
- Users should overload this method should their factor not default to `.Z<:ParametricType`.
- First design choice was to restrict this function to returning coordinates
  - See https://github.com/JuliaRobotics/RoME.jl/issues/465
  - Pay attention to which tangent space point is used for converting points on a manifold to coordinates,
    - Originally written just for Lie Groups to support legacy, but future needs may well alter the design.
- Original design driven by parametric solve and dead reckon tethering.

See also: [`accumulateFactorMeans`](@ref), [`solveFactorParametric`](@ref)
"""
function getMeasurementParametric end

function getMeasurementParametric(Z)
  return error(
    "$(typeof(Z)) is not supported, please use non-parametric or open an issue if it should be",
  )
end

function getMeasurementParametric(Z::Normal)
  meas = mean(Z)
  iσ = 1 / std(Z)^2
  return [meas], reshape([iσ], 1, 1)
end

function getMeasurementParametric(Z::MvNormal)
  meas = mean(Z)
  iΣ = invcov(Z)
  return meas, iΣ
end

# the point `p` on the manifold is the mean
function getMeasurementParametric(s::ManifoldPrior)
  meas = s.p
  iΣ = invcov(s.Z)
  return meas, iΣ
end

function getMeasurementParametric(s::AbstractFactor)
  if hasfield(typeof(s), :Z)
    Z = s.Z
  else
    error(
      "getMeasurementParametric(::$(typeof(s))) not defined, please add it, or use non-parametric, or open an issue for help.",
    )
  end

  return getMeasurementParametric(Z)
end

getMeasurementParametric(fct::DFGFactor) = getMeasurementParametric(getFactorType(fct))
getMeasurementParametric(dfg::AbstractDFG, flb::Symbol) = getMeasurementParametric(getFactor(dfg, flb))

# maybe rename getMeasurementParametric to something like getNormalDistributionParams or getMeanCov

# default to point on manifold
function getFactorMeasurementParametric(fac::AbstractPrior)
  M = getManifold(fac)
  ϵ = getPointIdentity(M)
  dims = manifold_dimension(M)
  Xc, iΣ = getMeasurementParametric(fac)
  X = get_vector(M, ϵ, Xc, DefaultOrthogonalBasis())
  meas = convert(typeof(ϵ), exp(M, ϵ, X))
  iΣ = convert(SMatrix{dims, dims}, iΣ)
  meas, iΣ
end
# default to point on tangent vector
function getFactorMeasurementParametric(fac::AbstractRelative)
  M = getManifold(fac)
  ϵ = getPointIdentity(M)
  dims = manifold_dimension(M)
  Xc, iΣ = getMeasurementParametric(fac)
  measX = convert(typeof(ϵ), get_vector(M, ϵ, Xc, DefaultOrthogonalBasis()))
  iΣ = convert(SMatrix{dims, dims}, iΣ)
  measX, iΣ
end

getFactorMeasurementParametric(fct::DFGFactor) = getFactorMeasurementParametric(getFactorType(fct))
getFactorMeasurementParametric(dfg::AbstractDFG, flb::Symbol) = getFactorMeasurementParametric(getFactor(dfg, flb))

## ================================================================================================
## Parametric solve with Mahalanobis distance - CalcFactor
## ================================================================================================

#TODO maybe remove with Mixture rework see #1504
getFactorMechanics(f::AbstractFactor) = f
getFactorMechanics(f::Mixture) = f.mechanics

function CalcFactorMahalanobis(fg, fct::DFGFactor)
  fac_func = getFactorType(fct)
  varOrder = getVariableOrder(fct)

  # NOTE, use getMeasurementParametric on DFGFactor{<:CCW} to allow special cases like OAS factors
  _meas, _iΣ = getFactorMeasurementParametric(fct) # fac_func
  
  # make sure its a tuple TODO Fix with mixture rework #1504
  meas = typeof(_meas) <: Tuple ? _meas : (_meas,)
  iΣ = typeof(_iΣ) <: Tuple ? _iΣ : (_iΣ,)

  cache = preambleCache(fg, getVariable.(fg, varOrder), getFactorType(fct))

  multihypo = getSolverData(fct).multihypo
  nullhypo = getSolverData(fct).nullhypo

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

  return CalcFactorMahalanobis(fct.label, getFactorMechanics(fac_func), cache, varOrder, meas, iΣ, special)
end

# This is where the actual parametric calculation happens, CalcFactor equivalent for parametric
# function (cfp::CalcFactorMahalanobis{FT, 1, C, MEAS, D, L, Nothing})(variables...) where {FT, C, MEAS, D, L, Nothing}# AbstractArray{T} where T <: Real
#   # call the user function
#   res = cfp.calcfactor!(cfp.meas..., variables...)
#   # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
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
    getFactorType(fct) isa MetaPrior ? continue : nothing
    calcFactors[fct.label] = CalcFactorMahalanobis(fg, fct)
  end
  return calcFactors
end

function getFactorTypesCount(facs::Vector{<:DFGFactor})
  typedict = OrderedDict{DataType, Int}()
  alltypes = OrderedDict{DataType, Vector{Symbol}}()
  for f in facs
    facType = typeof(getFactorType(f))
    cnt = get!(typedict, facType, 0)
    typedict[facType] = cnt + 1

    dt = get!(alltypes, facType, Symbol[])
    push!(dt, f.label)
  end
  #TODO tuple or vector?
  # vartypes = tuple(keys(typedict)...)
  factypes::Vector{DataType} = collect(keys(typedict))
  return factypes, typedict, alltypes
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

## ================================================================================================
## ================================================================================================
## New Parametric refactor WIP
## ================================================================================================
## ================================================================================================

## ================================================================================================
## LazyCase based on LazyBufferCache from PreallocationTools.jl
## ================================================================================================

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

## ================================================================================================
## GraphSolveStructures
## ================================================================================================

getVariableTypesCount(fg::AbstractDFG) = getVariableTypesCount(getVariables(fg))

function getVariableTypesCount(vars::Vector{<:DFGVariable})
  typedict = OrderedDict{DataType, Int}()
  alltypes = OrderedDict{DataType, Vector{Symbol}}()
  for v in vars
    varType = typeof(getVariableType(v))
    cnt = get!(typedict, varType, 0)
    typedict[varType] = cnt + 1

    dt = get!(alltypes, varType, Symbol[])
    push!(dt, v.label)
  end
  #TODO tuple or vector?
  # vartypes = tuple(keys(typedict)...)
  vartypes::Vector{DataType} = collect(keys(typedict))
  return vartypes, typedict, alltypes
end

buildGraphSolveManifold(fg::AbstractDFG) = buildGraphSolveManifold(getVariables(fg))

function buildGraphSolveManifold(vars::Vector{<:DFGVariable})
  vartypes, vartypecount, vartypeslist = getVariableTypesCount(vars)

  PMs = map(vartypes) do vartype
    N = vartypecount[vartype]
    G = getManifold(vartype)
    return NPowerManifold(G, N)
    # PowerManifold(G, NestedReplacingPowerRepresentation(), N)
    # PowerManifold(G, NestedPowerRepresentation(), N) #TODO investigate as it does not converge
  end
  M = ProductManifold(PMs...)
  return M, vartypes, vartypeslist
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
  Xc = get_coordinates(M, ϵ, X, DefaultOrthogonalBasis())
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
  # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
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

  # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)

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
      p[gsc.M, i][j] = getVariableSolverData(fg, vId, solveKey).val[1]
    end
  end
end

function _get_dim_ranges(dims::NTuple{N,Any}) where {N}
  dims_acc = accumulate(+, vcat(1, SVector(dims)))
  return ntuple(i -> (dims_acc[i]:(dims_acc[i] + dims[i] - 1)), Val(N))
end

#NOTE this only works with a product of power manifolds
function getComponentsCovar(@nospecialize(PM::ProductManifold), Σ::AbstractMatrix)
  dims = manifold_dimension.(PM.manifolds)
  dim_ranges = _get_dim_ranges(dims)

  subsigmas = map(zip(dim_ranges, PM.manifolds)) do v
    r = v[1]
    M = v[2]
    return _getComponentsCovar(M, view(Σ, r, r))
  end

  return ArrayPartition(subsigmas...)
end

function _getComponentsCovar(@nospecialize(PM::PowerManifold), Σ::AbstractMatrix)
  M = PM.manifold
  dim = manifold_dimension(M)
  subsigmas = map(Manifolds.get_iterator(PM)) do i
    r = ((i - 1) * dim + 1):(i * dim)
    return Σ[r, r]
  end

  return subsigmas
end

function _getComponentsCovar(@nospecialize(PM::NPowerManifold), Σ::AbstractMatrix)
  M = PM.manifold
  dim = manifold_dimension(M)
  subsigmas = map(Manifolds.get_iterator(PM)) do i
    r = ((i - 1) * dim + 1):(i * dim)
    return Σ[r, r]
  end

  return subsigmas
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
  get_coordinates!(M, Xc, ϵ, X, DefaultOrthogonalBasis())

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

## Original
# ==============================

function _totalCost(fg, cfdict::OrderedDict{Symbol, <:CalcFactorMahalanobis}, flatvar, Xc)
  #
  obj = zero(eltype(Xc))
  for (fid, cfp) in cfdict
    varOrder = cfp.varOrder

    Xparams = [
      getPoint(getVariableType(fg, varId), view(Xc, flatvar.idx[varId])) for
      varId in varOrder
    ]

    # call the user function
    # retval = cfp(Xparams...)
    res = cfp(cfp.meas..., Xparams...)
    # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
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
  varIds = [frontals; separators]

  sfg = issetequal(varIds, listVariables(fg)) ? fg : buildSubgraph(fg, varIds, 1)

  flatvar = FlatVariables(fg, varIds)

  for vId in varIds
    p = getVariableSolverData(fg, vId, solvekey).val[1]
    flatvar[vId] = getCoordinates(getVariableType(fg, vId), p)
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
    p = getPoint(getVariableType(fg, key), rv[r])
    push!(d, key => (val = p, cov = Σ[r, r]))
  end

  return (opti = d, stat = result, varIds = flatvar.idx, Σ = Σ)
end

## ================================================================================================
## UNDER DEVELOPMENT Parametric solveTree utils
## ================================================================================================

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

"""
    $SIGNATURES

"""
function calculateCoBeliefMessage(soldict, Σ, flatvars, separators, frontals)
  Aidx = IIF.collectIdx(flatvars, separators)
  Cidx = IIF.collectIdx(flatvars, frontals)

  #marginalize separators
  A = Σ[Aidx, Aidx]
  #marginalize frontals
  C = Σ[Cidx, Cidx]
  # cross
  B = Σ[Aidx, Cidx]

  Σₘ = deepcopy(A)
  if length(separators) == 0
    return (varlbl = Symbol[], μ = Float64[], Σ = Matrix{Float64}(undef, 0, 0))

  elseif length(separators) == 1

    # create messages
    return (varlbl = deepcopy(separators), μ = soldict[separators[1]].val, Σ = A)

  elseif length(separators) == 2
    A = Σₘ[1, 1]
    C = Σₘ[2, 2]
    B = Σₘ[1, 2]

    #calculate covariance between separators
    ΣA_B = A - B * inv(C) * B'
    # create messages
    m2lbl = deepcopy(separators)
    m2cov = isa(ΣA_B, Matrix) ? ΣA_B : fill(ΣA_B, 1, 1)
    m2val = soldict[m2lbl[2]].val - soldict[m2lbl[1]].val
    return (varlbl = m2lbl, μ = m2val, Σ = m2cov)

  else
    error("Messages with more than 2 seperators are not supported yet")
  end
end

## ================================================================================================
## Parametric utils
## ================================================================================================

## SANDBOX of usefull development functions to be cleaned up
"""
    $SIGNATURES
Update the parametric solver data value and covariance.
"""
function updateSolverDataParametric! end

function updateSolverDataParametric!(
  vnd::VariableNodeData,
  val::AbstractArray,
  cov::AbstractMatrix,
)
  # fill in the variable node data value
  vnd.val[1] = val
  #calculate and fill in covariance
  vnd.bw .= cov
  return vnd
end

function updateSolverDataParametric!(
  v::DFGVariable,
  val::AbstractArray,
  cov::AbstractMatrix;
  solveKey::Symbol = :parametric,
)
  vnd = getSolverData(v, solveKey)
  return updateSolverDataParametric!(vnd, val, cov)
end


"""
    $SIGNATURES
Add parametric solver to fg, batch solve using [`solveGraphParametric`](@ref) and update fg.
"""
function solveGraphParametricOptim!(
  fg::AbstractDFG; 
  init::Bool = true, 
  solveKey::Symbol = :parametric, # FIXME, moot since only :parametric used for parametric solves
  initSolveKey::Symbol = :default, 
  verbose = false,
  kwargs...
)
  # make sure variables has solverData, see #1637
  makeSolverData!(fg; solveKey)
  if !(:parametric in fg.solverParams.algorithms)
    addParametricSolver!(fg; init = init)
  elseif init
    initParametricFrom!(fg, initSolveKey; parkey=solveKey)
  end

  vardict, result, varIds, Σ = solveGraphParametricOptim(fg; verbose, kwargs...)

  updateParametricSolution!(fg, vardict)

  return vardict, result, varIds, Σ
end

"""
    $SIGNATURES
Initialize the parametric solver data from a different solution in `fromkey`.

DevNotes
- TODO, keyword `force` not wired up yet.
"""
function initParametricFrom!(
  fg::AbstractDFG,
  fromkey::Symbol = :default;
  parkey::Symbol = :parametric,
  onepoint = false,
  force::Bool = false,
)
  #
  if onepoint
    for v in getVariables(fg)
      fromvnd = getSolverData(v, fromkey)
      dims = getDimension(v)
      getSolverData(v, parkey).val[1] = fromvnd.val[1]
      getSolverData(v, parkey).bw[1:dims, 1:dims] = LinearAlgebra.I(dims)
    end
  else
    for var in getVariables(fg)
      dims = getDimension(var)
      μ, Σ = calcMeanCovar(var, fromkey)
      getSolverData(var, parkey).val[1] = μ
      getSolverData(var, parkey).bw[1:dims, 1:dims] = Σ
    end
  end
end

"""
    $SIGNATURES
Add the parametric solveKey to all the variables in fg if it doesn't exists.
"""
function addParametricSolver!(fg; init = true)
  if !(:parametric in fg.solverParams.algorithms)
    push!(fg.solverParams.algorithms, :parametric)
    foreach(
      v -> IIF.setDefaultNodeDataParametric!(v, getVariableType(v); initialized = false),
      getVariables(fg),
    )
    if init
      initParametricFrom!(fg)
    end
  else
    error("parametric solvekey already exists")
  end
  return nothing
end

"""
    $SIGNATURES
Update the fg from solution in vardict and add MeanMaxPPE (all just mean). Usefull for plotting
"""
function updateParametricSolution!(sfg, vardict::AbstractDict; solveKey::Symbol = :parametric)
  for (v, val) in vardict
    vnd = getSolverData(getVariable(sfg, v), solveKey)
    # Update the variable node data value and covariance
    updateSolverDataParametric!(vnd, val.val, val.cov)
    #fill in ppe as mean
    Xc = collect(getCoordinates(getVariableType(sfg, v), val.val))
    ppe = MeanMaxPPE(solveKey, Xc, Xc, Xc)
    getPPEDict(getVariable(sfg, v))[solveKey] = ppe
  end
end

function updateParametricSolution!(fg, M, labels::AbstractArray{Symbol}, vals, Σ; solveKey::Symbol = :parametric)
  
  if !isnothing(Σ)
    covars = getComponentsCovar(M, Σ)
  end

  for (i, (v, val)) in enumerate(zip(labels, vals))
    vnd = getSolverData(getVariable(fg, v), solveKey)
    covar = isnothing(Σ) ? vnd.bw : covars[i]
    # Update the variable node data value and covariance
    updateSolverDataParametric!(vnd, val, covar)#FIXME add cov
    #fill in ppe as mean
    Xc = collect(getCoordinates(getVariableType(fg, v), val))
    ppe = MeanMaxPPE(solveKey, Xc, Xc, Xc)
    getPPEDict(getVariable(fg, v))[solveKey] = ppe
  end

end

function createMvNormal(val, cov)
  #TODO do something better for properly formed covariance, but for now just a hack...FIXME
  if all(diag(cov) .> 0.001) && isapprox(cov, transpose(cov); rtol = 1e-4)
    return MvNormal(val, Symmetric(cov))
  else
    @error("Covariance matrix error", cov)
    # return nothing # FIXME, blanking nothing during #459 consolidation
    return MvNormal(val, ones(length(val)))
  end
end

function createMvNormal(v::DFGVariable, key = :parametric)
  if key == :parametric
    vnd = getSolverData(v, :parametric)
    dims = vnd.dims
    return createMvNormal(vnd.val[1:dims, 1], vnd.bw[1:dims, 1:dims])
  else
    @warn "Trying MvNormal Fit, replace with PPE fits in future"
    return fit(MvNormal, getSolverData(v, key).val)
  end
end

#TODO this is still experimental and a POC
function getInitOrderParametric(fg; startIdx::Symbol = lsfPriors(fg)[1])
  order = DFG.traverseGraphTopologicalSort(fg, startIdx)
  filter!(order) do l
    return isVariable(fg, l)
  end
  return order
end

function autoinitParametricOptim!(
  fg,
  varorderIds = getInitOrderParametric(fg);
  reinit = false,
  algorithm = Optim.NelderMead,
  algorithmkwargs = (initial_simplex = Optim.AffineSimplexer(0.025, 0.1),),
  kwargs...
)
  @showprogress for vIdx in varorderIds
    autoinitParametricOptim!(fg, vIdx; reinit, algorithm, algorithmkwargs, kwargs...)
  end
  return nothing
end

function autoinitParametricOptim!(dfg::AbstractDFG, initme::Symbol; kwargs...)
  return autoinitParametricOptim!(dfg, getVariable(dfg, initme); kwargs...)
end

function autoinitParametricOptim!(
  dfg::AbstractDFG,
  xi::DFGVariable;
  solveKey = :parametric,
  reinit::Bool = false,
  kwargs...,
)
  #

  initme = getLabel(xi)
  vnd = getSolverData(xi, solveKey)
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
    #fill in ppe as mean
    Xc = collect(getCoordinates(getVariableType(xi), val))
    ppe = MeanMaxPPE(:parametric, Xc, Xc, Xc)
    getPPEDict(xi)[:parametric] = ppe

    # updateVariableSolverData!(dfg, xi, solveKey, true; warn_if_absent=false)    
    # updateVariableSolverData!(dfg, xi.label, getSolverData(xi, solveKey), :graphinit, true, Symbol[]; warn_if_absent=false)
  else
    result = nothing
  end

  return result#isInitialized(xi, solveKey)
end

#
