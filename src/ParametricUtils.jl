## ================================================================================================
## ArrayPartition getPointIdentity (identity_element)
## ================================================================================================
import Manifolds: TraitList, IsExplicitDecorator, EmptyTrait, AdditionGroupTrait, trait, _padpoint!, next_trait
using RecursiveArrayTools: ArrayPartition

import DistributedFactorGraphs: getPointIdentity


function getPointIdentity(G::ProductGroup, numtype=Float64)
  M = G.manifold
  return ProductRepr(map(x->getPointIdentity(x, numtype), M.manifolds))
  # return ArrayPartition(map(x->getPointIdentity(x, numtype), M.manifolds))
end

function getPointIdentity(G::SpecialOrthogonal{N}, numtype=Float64) where N
  return SMatrix{N,N, numtype}(I)
  # return Matrix(one(numtype)*I, N, N)
end

function getPointIdentity(G::TranslationGroup{Tuple{N}}, numtype=Float64) where N
  return zeros(SVector{N,numtype})
  # return zeros(numtype, N)
end

# fallback 
function getPointIdentity(G::GroupManifold, numtype=Float64)
  return error("not implemented")
end

function getPointIdentity(G::ProductManifold, numtype=Float64)
  return ProductRepr(map(x->getPointIdentity(x,numtype), G.manifolds))
end

function getPointIdentity(M::Manifolds.PowerManifoldNestedReplacing, numtype=Float64)
  N = Manifolds.get_iterator(M).stop
  return fill(getPointIdentity(M.manifold, numtype), N)
end

function getPointIdentity(M::PowerManifold, numtype=Float64)
  N = Manifolds.get_iterator(M).stop
  return fill(getPointIdentity(M.manifold, numtype), N)
end

function getPointIdentity(G::SemidirectProductGroup, numtype=Float64)
  M = base_manifold(G)
  N, H = M.manifolds
  np = getPointIdentity(N,numtype)
  hp = getPointIdentity(H,numtype)
  # return ArrayPartition(np, hp)
  return ProductRepr(np, hp)
end

function getPointIdentity(::SpecialOrthogonal{2})
  return SA[1.0 0.0; 0.0 1.0]
end

function getPointIdentity(::SpecialOrthogonal{3})
  return SA[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
end

function getPointIdentity(::SpecialEuclidean{2})
  return ProductRepr(SA[0.,0.], SA[1.0 0.0; 0.0 1.0])
  # return ArrayPartition(SA[0.,0.], SA[1.0 0.0; 0.0 1.0])
end

function getPointIdentity(::SpecialEuclidean{3})
  # or ArrayPartition
  return ProductRepr(SA[0.,0.,0.], SA[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
end


function Manifolds.allocate_result(G::SemidirectProductGroup, ::typeof(getPointIdentity))
  M = base_manifold(G)
  N, H = M.manifolds
  np = allocate_result(N, getPointIdentity)
  hp = allocate_result(H, getPointIdentity)
  # return ArrayPartition(np, hp)
  return ProductRepr(np, hp)
end

##

# test
# M = ProductManifold(SpecialEuclidean(2), SpecialEuclidean(3), SpecialOrthogonal(3), TranslationGroup(2));
# G = ProductGroup(M);
# getPointIdentity(G)

## ================================================================================================
## GraphSolveStructures
## ================================================================================================
using OrderedCollections

function getVariableTypesCount(fg::AbstractDFG)

  vars = getVariables(fg)
  typedict = OrderedDict{DataType, Int}()
  alltypes = OrderedDict{DataType,Vector{Symbol}}()
  for v in vars
      varType = typeof(getVariableType(v))
      cnt = get!(typedict, varType, 0)
      typedict[varType] = cnt+1

      dt = get!(alltypes, varType, Symbol[])
      push!(dt, v.label)
  end
  vartypes = tuple(keys(typedict)...)
  return vartypes, typedict, alltypes
end

function buildGraphSolveManifold(fg)
  vartypes, vartypecount, vartypeslist = getVariableTypesCount(fg)
  M = mapreduce(ProductManifold, vartypes) do vartype
      N = vartypecount[vartype] 
      G = getManifold(vartype)
      PowerManifold(G, NestedReplacingPowerRepresentation(), N)
      # PowerManifold(G, NestedPowerRepresentation(), N)
  end
  return M, vartypes, vartypeslist
end


# function Manifolds.vee(M::AbstractPowerManifold, X)
#   rep_size = representation_size(M.manifold)
#   p = getPointIdentity(M.manifold, Float64)
#   vs = [
#       vee(M.manifold, p, Manifolds._read(M, rep_size, X, i)) for i in Manifolds.get_iterator(M)
#   ]
#   return reduce(vcat, reshape(vs, length(vs)))
# end

# function Manifolds.vee(M::ProductManifold, X)
#     reps = map(
#         vee,
#         M.manifolds,
#         submanifold_components(M, X))
#     return vcat(reps...)
# end

# function Manifolds.log!(M::ProductManifold, X, ::Identity{ProductOperation}, q)
#   map(
#       log!,
#       M.manifolds,
#       submanifold_components(M, X),
#       submanifold_components(M, p),
#       submanifold_components(M, q),
#   )
#   return X
# end

# function Manifolds.hat!(M::Manifolds.PowerManifoldNestedReplacing, Y, p, c)
#   B=VeeOrthogonalBasis()
#   dim = manifold_dimension(M.manifold)
#   rep_size = representation_size(M.manifold)
#   v_iter = 1
#   for i in Manifolds.get_iterator(M)
#      get_vector!(
#           M.manifold,
#           Manifolds._write(M, rep_size, Y, i),
#           Manifolds._read(M, rep_size, Y, i),
#           c[v_iter:(v_iter + dim - 1)],
#           B,
#       )
#       v_iter += dim
#   end
#   return Y
# end


## ================================================================================================
## FlatVariables - used for packing variables for optimization
## ================================================================================================

struct FlatVariables{T<:Real}
  X::Vector{T}
  p::ProductRepr
  M::ProductGroup
  idx::Dict{Symbol, UnitRange{Int}}
end

function FlatVariables(fg::AbstractDFG, varIds::Vector{Symbol})

  index = 1
  idx = Dict{Symbol, UnitRange{Int}}()
  for vid = varIds
    v = getVariable(fg, vid)
    dims = getDimension(v)
    idx[vid] = index:(index+dims-1)
    index += dims
  end
  M = ProductGroup(ProductManifold(getManifold.(fg, varIds)...))
  p = ProductRepr(getPointIdentity.(getVariableType.(fg,varIds))...)
  return FlatVariables(Vector{Float64}(undef, index-1), p, M, idx)
end

function Base.setindex!(flatVar::FlatVariables{T}, val::Vector{T}, vId::Symbol) where T<:Real
  if length(val) == length(flatVar.idx[vId])
    flatVar.X[flatVar.idx[vId]] .= val
  else
    error("array could not be broadcast to match destination")
  end
end

function Base.getindex(flatVar::FlatVariables{T}, vId::Symbol) where T<:Real
  return flatVar.X[flatVar.idx[vId]]
end


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
LazyCache(f::F = allocate) where {F <: Function} = LazyCache(Dict{Tuple{DataType, Symbol}, Any}(), f)

# override the [] method
function Base.getindex(cache::LazyCache, u::T, varname::Symbol) where T
    val = get!(cache.dict, (T, varname)) do 
      cache.fnc(u)
    end::T
    return val
end



function getCoordCache!(cache::LazyCache, M, T::DataType, varname::Symbol)
  val = get!(cache.dict, (T, varname)) do 
    Vector{T}(undef,manifold_dimension(M))
  end::Vector{T}
  return val
end



## ================================================================================================
## Parametric Factors
## ================================================================================================

"""
    $SIGNATURES

Returns the parametric measurement for a factor as a tuple (measurement, inverse covariance) for parametric inference (assuming Gaussian).
Defaults to find the parametric measurement at field `Z`, fields `Zij` and `z` are deprecated for standardization.

Notes
- Users should overload this method should their factor not default to `.Z<:ParametricType`.
- First design choice was to restrict this function to returning coordinates
  - See https://github.com/JuliaRobotics/RoME.jl/issues/465
  - Pay attention to which tangent space point is used for converting points on a manifold to coordinates,
    - Originally written just for Lie Groups to support legacy, but future needs may well alter the design.
- Original design driven by parametric solve and dead reckon tethering.

See also: [`accumulateFactorMeans`](@ref), [`solveFactorParameteric`](@ref)
"""
function getMeasurementParametric end

function getMeasurementParametric(Z)
  error("$(typeof(Z)) is not supported, please use non-parametric or open an issue if it should be")
end

function getMeasurementParametric(Z::Normal)
  meas = mean(Z)
  iσ = 1/std(Z)^2
  return [meas], reshape([iσ],1,1)
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
    @info "getMeasurementParametric falls back to using field `.Z` by default. Extend it for more complex factors." maxlog=1
  else
    error("getMeasurementParametric(::$(typeof(s))) not defined, please add it, or use non-parametric, or open an issue for help.")
  end
  
  return getMeasurementParametric(Z)
end



## ================================================================================================
## Parametric solve with Mahalanobis distance - CalcFactor
## ================================================================================================



getFactorMechanics(f::AbstractFactor) = f
getFactorMechanics(f::Mixture) = f.mechanics

function CalcFactorMahalanobis(fg, fct::DFGFactor)
  cf = getFactorType(fct)
  varOrder = getVariableOrder(fct)
  varTypes = getVariableType.(fg, varOrder) .|> typeof
  _meas, _iΣ = getMeasurementParametric(cf)
  M = getManifold(getFactorType(fct))

  (typeof(_meas) <: Tuple) && error("Not implemented $_meas")

  ϵ = getPointIdentity(M)
  if cf isa AbstractPrior
    meas = (exp(M, ϵ, hat(M, Identity(M), _meas)),)
  else  
    meas = (hat(M, Identity(M), _meas),)
  end

  iΣ = typeof(_iΣ) <: Tuple ? _iΣ : (_iΣ,)

  cache = preambleCache(fg, getVariable.(fg, varOrder), getFactorType(fct))
  cache = isnothing(cache) ? NamedTuple() : cache

  calcf = CalcFactor(getFactorMechanics(cf), nothing, 0, 0, nothing, nothing, true, nothing)
  # calcf = CalcFactor(getFactorMechanics(cf), _getFMdThread(fct), 0, 0, nothing, nothing, true, cache)

  multihypo = getSolverData(fct).multihypo
  nullhypo = getSolverData(fct).nullhypo

  # FIXME, type instability, use dispatch instead of if-else
  if length(multihypo) > 0
    special = MaxMultihypo(multihypo)
  elseif nullhypo > 0 
    special = MaxNullhypo(nullhypo)
  elseif cf isa Mixture
    special = MaxMixture(cf.diversity.p, Ref(0))
  else
    special = nothing
  end

  return CalcFactorMahalanobis(calcf, varOrder, varTypes, meas, iΣ, special)
end

# This is where the actual parametric calculation happens, CalcFactor equivalent for parametric
function (cfp::CalcFactorMahalanobis)(variables...)
  # call the user function (be careful to call the new CalcFactor version only!!!)
  # @warn "1" a=typeof(cfp.meas[1])

  res = cfp.calcfactor!(cfp.meas[1], variables...)
  # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
  return res' * cfp.iΣ[1] * res
end

# function (cfp::CalcFactorMahalanobis{CalcFactor{T, U, V, W, C}})(variables...) where {T<:AbstractFactor, U, V, W, C}
#   # call the user function (be careful to call the new CalcFactor version only!!!)
#   res = cfp.calcfactor!(cfp.meas[1], variables...)
#   # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
#   return res' * cfp.iΣ[1] * res
# end

# function (cfp::CalcFactorMahalanobis{CalcFactor{T, U, V, W, C}})(variables...) where {T<:AbstractPrior, U, V, W, C}
#   # TODO should prior functions follow the same factor definition with a measurement on tangant?
#   M = getManifold(cfp.calcfactor!.factor)

#   ϵ = getPointIdentity(M)
#   m = exp(M, ϵ, cfp.meas[1])  

#   res = cfp.calcfactor!(m, variables...)
#   # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
#   return res' * cfp.iΣ[1] * res
# end

# function (cfp::CalcFactorMahalanobis{CalcFactor{T, U, V, W, C}})(variables...) where {T<:Union{ManifoldFactor, ManifoldPrior} , U, V, W, C}
#   # call the user function (be careful to call the new CalcFactor version only!!!)
#   M = cfp.calcfactor!.factor.M
#   X = cfp.calcfactor!(cfp.meas[1], variables...)
#   # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
#   # return mahalanobus_distance2(M, X, cfp.iΣ[1])
  
#   #TODO do something about basis?
#   # Xc = get_coordinates(M, variables[1], X, DefaultOrthogonalBasis())
#   # Xc = get_coordinates(M, variables[1], X, DefaultOrthonormalBasis())
#   Xc = vee(M, variables[1], X)
#   return X, Xc' * cfp.iΣ[1] * Xc
# end

function calcFactorMahalanobisDict(fg)
  calcFactors = OrderedDict{Symbol, CalcFactorMahalanobis}()
  for fct in getFactors(fg)
    calcFactors[fct.label] = CalcFactorMahalanobis(fg, fct)
  end
  return calcFactors
end


cost_cfp(cfp::CalcFactorMahalanobis, M, p, varOrder::NTuple{1,Int}) = cfp(p[M, varOrder[1]])
cost_cfp(cfp::CalcFactorMahalanobis, M, p, varOrder::NTuple{2,Int}) = cfp(p[M, varOrder[1]], p[M, varOrder[2]])
cost_cfp(cfp::CalcFactorMahalanobis, M, p, varOrder::NTuple{3,Int}) = cfp(p[M, varOrder[1]], p[M, varOrder[2]], p[M, varOrder[3]])

cost_cfp(cfp::CalcFactorMahalanobis, M, p, vi::NTuple{1,Tuple{Int64, Int64}}) = cfp(p[M, vi[1][1]][vi[1][2]])
cost_cfp(cfp::CalcFactorMahalanobis, M, p, vi::NTuple{2,Tuple{Int64, Int64}}) = cfp(p[M, vi[1][1]][vi[1][2]], p[M, vi[2][1]][vi[2][2]])

function _totalCost(fg, 
                    cfdict::OrderedDict{Symbol, <:CalcFactorMahalanobis},
                    flatvar,
                    point_varOrders,
                    Xc )
  #
  
  
  if (test_new=true)
    M = flatvar.M
    # p = flatvar.p
    # hat!(M, p, Identity(M), Xc)
    # exp!(M, p, identity_element(M), p)

    p = exp(M, identity_element(M), hat(M, Identity(M), Xc))

  end

  
  #TODO create once fur multithread
  # obj = zero(eltype(Xc))
  obj = zeros(eltype(Xc),(Threads.nthreads()))

  Threads.@threads for fid in collect(keys(cfdict))
    cfp = cfdict[fid]
  # for (fid, cfp) in cfdict 
    # call the user function
    if test_new
      retval = cost_cfp(cfp, M, p, point_varOrders[fid])
    else
      varOrder = cfp.varOrder
      Xparams = [getPoint(getVariableType(fg, varId), view(Xc, flatvar.idx[varId])) for varId in varOrder]
      retval = cfp(Xparams...)
    end    
    
    obj[Threads.threadid()] += retval
    # obj += retval
  end

  return 1/2*sum(obj)
  
  # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
  # return 1/2*obj
end

## :forward test_new==false
# 57.366395 seconds (33.10 M allocations: 1.626 GiB, 0.65% gc time, 99.43% compilation time)
# 0.309763 seconds (1.42 M allocations: 321.207 MiB, 16.38% gc time)

## :forward test_new==true N=30
# 5.262519 seconds (15.64 M allocations: 963.398 MiB, 4.91% gc time, 94.58% compilation time)
# 0.230315 seconds (758.36 k allocations: 244.261 MiB, 9.56% gc time)

## :forward test_new==true N=50
# 7.334445 seconds (13.06 M allocations: 1.619 GiB, 4.52% gc time, 80.12% compilation time)
# 1.215492 seconds (3.77 M allocations: 1.138 GiB, 11.52% gc time)

# 1 thread
# 1.154163 seconds (3.77 M allocations: 488.173 MiB, 12.62% gc time)
# 1.109158 seconds (3.77 M allocations: 488.173 MiB, 8.16% gc time)
# 4 threads
# 0.727205 seconds (3.86 M allocations: 498.374 MiB, 13.44% gc time)
# 0.739042 seconds (3.86 M allocations: 498.372 MiB, 17.42% gc time)

# function cacheByType!(cache::LazyCache, T::DataType, varname::Symbol, fnc, x...) 
#   val = get!(cache.dict, (T, varname)) do 
#     fnc(x...)
#   end
#   return val
# end

# if numtype in cashedTypes
#   ϵ = cacheByType!(lazycache, numtype, :ϵ, identity, 0)
#   p = cacheByType!(lazycache, numtype, :p, identity, 0)
#   X = cacheByType!(lazycache, numtype, :X, identity, 0)
# else
#   @warn "cache miss"
#   ϵ = cacheByType!(lazycache, numtype, :ϵ, getPointIdentity, M, numtype)
#   p = cacheByType!(lazycache, numtype, :p, deepcopy, ϵ)
#   X = cacheByType!(lazycache, numtype, :X, deepcopy, ϵ)
#   @warn "added to cache"
# end

# @profview foreach((_Xc)->IIF._totalCost2(fg, cfd, gsc, varOrderDict, _Xc), [Xc+randn(length(Xc)) for _=1:10])


struct GraphSolveBuffers{T<:Real, U}
  ϵ::U
  p::U
  X::U
  Xc::Vector{T}
end

function GraphSolveBuffers(M, ::Type{T}) where T
  ϵ = getPointIdentity(M, T)
  p = deepcopy(ϵ)# allocate_result(M, getPointIdentity)
  X = deepcopy(ϵ) #allcoate(p)
  Xc = get_coordinates(M, ϵ, X, DefaultOrthogonalBasis())
  GraphSolveBuffers(ϵ, p, X, Xc)
end

struct GraphSolveContainer
  M::AbstractManifold # ProductManifold or ProductGroup
  buffers::OrderedDict{DataType, GraphSolveBuffers}
  varTypes::Vector{DataType}
  varTypesIds::OrderedDict{DataType, Vector{Symbol}}
  cfdict::OrderedDict{Symbol, CalcFactorMahalanobis}
  varOrderDict::OrderedDict{Symbol, Tuple{Tuple{Int,Int}, Vararg{Tuple{Int,Int}}}}
end


function GraphSolveContainer(fg)

  M, varTypes, varTypesIds = buildGraphSolveManifold(fg)
  buffs = OrderedDict{DataType,GraphSolveBuffers}()
  cfd = calcFactorMahalanobisDict(fg)

  varOrderDict = OrderedDict{Symbol, Tuple{Tuple{Int,Int}, Vararg{Tuple{Int,Int}}}}()#NTuple{N,Tuple{Int64, Int64}}}()
  #TODO find better way for indexing
  for (fid, cfp) in cfd 
    varOrder = cfp.varOrder
    prod_pow_idx = map(varOrder) do v
      v_type = getVariableType(fg, v) |> typeof
      prod_idx = findfirst(==(v_type), varTypes)
      vars = varTypesIds[varTypes[prod_idx]]
      pow_idx = findfirst(==(v), vars)
      (prod_idx, pow_idx)
    end
    prod_pow_idx = tuple(prod_pow_idx...)
    varOrderDict[fid] = prod_pow_idx
  end

  return GraphSolveContainer(M, buffs, varTypes, varTypesIds, cfd, varOrderDict)
end

function getGraphSolveCache!(gsc::GraphSolveContainer, ::Type{T}) where T<:Real
  cache = gsc.buffers
  M = gsc.M
  val = get!(cache, T) do 
    @info "cache miss, cacheing" T
    GraphSolveBuffers(M, T)
  end
  return val
end

function _toPoints2!(M::AbstractManifold, buffs::GraphSolveBuffers{T,U}, Xc::Vector{T}) where {T,U}
  ϵ = buffs.ϵ
  p = buffs.p
  X = buffs.X
  get_vector!(M, X, ϵ, Xc, DefaultOrthogonalBasis())
  exp!(M, p, ϵ, X)
  return p::U
end

Base.firstindex(::OrderedCollections.OrderedDict{Symbol, IncrementalInference.CalcFactorMahalanobis}) = error()

function (gsc::GraphSolveContainer)(Xc::Vector{T}) where T <: Real
  #
  buffs = getGraphSolveCache!(gsc, T)

  cfdict = gsc.cfdict
  varOrderDict = gsc.varOrderDict

  M = gsc.M 

  p = _toPoints2!(M, buffs, Xc)

  #NOTE multi threaded option
  obj = zeros(T,(Threads.nthreads()))
  Threads.@threads for fid in collect(keys(cfdict))
    cfp = cfdict[fid]

  #NOTE single thread option
  # obj::T = zero(T)
  # for (fid, cfp) in cfdict 

    prod_pow_idx = varOrderDict[fid]

    # call the user function
    retval = cost_cfp(cfp, M, p, prod_pow_idx)

    #NOTE multi threaded option
    obj[Threads.threadid()] += retval
    # NOTE single thread option
    # obj += retval
  end

  # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)

  #NOTE multi threaded option
  return sum(obj)/2
  # NOTE single thread option
  # return obj/2
end

#fg = generateCanonicalFG_Honeycomb!()

  # copy variables from graph
function initPoints!(p, gsc, fg::AbstractDFG, solveKey=:parametric)
  for (i, vartype) in enumerate(gsc.varTypes)
    varIds = gsc.varTypesIds[vartype]
    for (j,vId) in enumerate(varIds)
      p[gsc.M, i][j] = getVariableSolverData(fg, vId, solveKey).val[1]
    end
  end
end


# fg = generateGraph_Hexagonal(;fg=initfg(GraphsDFG;solverParams=SolverParams(algorithms=[:default, :parametric])), graphinit=false)
# 0.082069 seconds (517.90 k allocations: 56.727 MiB)
# (770.17 k allocations: 136.930 MiB, 11.95% gc time)
using ForwardDiff
function solveGraphParametric2(fg::AbstractDFG;
                              computeCovariance::Bool = false,
                              solveKey::Symbol=:parametric,
                              autodiff = :forward,
                              algorithm=Optim.BFGS,
                              algorithmkwargs=(), # add manifold to overwrite computed one
                              options = Optim.Options(allow_f_increases=true,
                                                      time_limit = 100,
                                                      # show_trace = true,
                                                      # show_every = 1,
                                                      ))
# 
  # Build the container  
  # fg = Main.generateGraph_Hexagonal(;fg=initfg(GraphsDFG;solverParams=SolverParams(algorithms=[:default, :parametric])), graphinit=false)                                                    
  # T = Optim.ForwardDiff.Dual{ForwardDiff.Tag{IncrementalInference.var"#tc#667"{GraphsDFG{SolverParams, DFGVariable, DFGFactor}, OrderedDict{Symbol, Tuple{Tuple{Int64, Int64}, Vararg{Tuple{Int64, Int64}}}}, OrderedDict{Symbol, IncrementalInference.CalcFactorMahalanobis}, IncrementalInference.LazyGraphSolveCache}, Float64}, Float64, 12}
  # lazyCache = LazyGraphSolveCache()
  gsc = GraphSolveContainer(fg)
  # buffs = getGraphSolveCache!(gsc, ForwardDiff.Dual{ForwardDiff.Tag{IncrementalInference.GraphSolveContainer, Float64}, Float64, 12})
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
  initValues .+= randn(length(Xc))*0.0001

  #optim setup and solve
  alg = algorithm(; algorithmkwargs...)
  tdtotalCost = Optim.TwiceDifferentiable(gsc, initValues, autodiff = autodiff)

  @time result = Optim.optimize(tdtotalCost, initValues, alg, options)
  rv = Optim.minimizer(result)

  # optionally compute hessian for covariance
  Σ = if computeCovariance
    H = Optim.hessian!(tdtotalCost, rv)
    pinv(H)
  else
    nothing
  end

  #TODO better return 

  # d = OrderedDict{Symbol,NamedTuple{(:val, :cov),Tuple{Vector{Float64},Matrix{Float64}}}}()

  # for key in varIds
  #   r = flatvar.idx[key]
  #   push!(d,key=>(val=rv[r],cov=Σ[r,r]))
  # end

  get_vector!(M, X, ϵ, rv, DefaultOrthogonalBasis())
  exp!(M, p, ϵ, X)

  # return d, result, flatvar.idx, Σ
  return (result=result, p=p, gsc=gsc, Σ=Σ)
end


export solveGraphParametric
"""
    $SIGNATURES

Batch solve a Gaussian factor graph using Optim.jl. Parameters can be passed directly to optim.
Notes:
  - Only :Euclid and :Circular manifolds are currently supported, own manifold are supported with `algorithmkwargs` (code may need updating though)
"""
function solveGraphParametric(fg::AbstractDFG;
                              computeCovariance::Bool = false,
                              solvekey::Symbol=:parametric,
                              autodiff = :forward,
                              algorithm=Optim.BFGS,
                              algorithmkwargs=(), # add manifold to overwrite computed one
                              options = Optim.Options(allow_f_increases=true,
                                                      time_limit = 100,
                                                      # show_trace = true,
                                                      # show_every = 1,
                                                      ))

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

  #TODO maybe remove sorting, just for convenience
  sort!(varIds, lt=natural_lt)

  flatvar = FlatVariables(fg, varIds)

  for vId in varIds
    p = getVariableSolverData(fg, vId, solvekey).val[1]
    flatvar[vId] = getCoordinates(getVariableType(fg,vId), p)
  end

  initValues = flatvar.X
  initValues .+= randn(length(initValues))*0.001

  alg = algorithm(; algorithmkwargs...)

  cfd = calcFactorMahalanobisDict(fg)
  
  
  point_varOrders = OrderedDict(fac_lbl => tuple(indexin(cfm.varOrder, varIds)...)  for (fac_lbl, cfm) in cfd)

  #TODO check if closure is correct for performance
  tc(x)=_totalCost(fg, cfd, flatvar, point_varOrders, x)

  # TODO remove, only testing
  # @time tc(initValues)

  tdtotalCost = Optim.TwiceDifferentiable(tc, initValues, autodiff = autodiff)

  @time result = Optim.optimize(tdtotalCost, initValues, alg, options)
  rv = Optim.minimizer(result)

  Σ = if computeCovariance
    H = Optim.hessian!(tdtotalCost, rv)
    pinv(H)
  else
    nothing
  end

  d = OrderedDict{Symbol,NamedTuple{(:val, :cov),Tuple{Vector{Float64},Matrix{Float64}}}}()

  for key in varIds
    r = flatvar.idx[key]
    cov = if Σ === nothing
      [0.;;]
    else
      Σ[r,r]
    end
    push!(d,key=>(val=rv[r],cov=cov))
  end

  return d, result, flatvar.idx, Σ
end

#TODO maybe consolidate with solveGraphParametric
#TODO WIP
```
    $SIGNATURES
Solve for frontal values only with values in seprarators fixed
```
function solveConditionalsParametric(fg::AbstractDFG,
                                    frontals::Vector{Symbol};
                                    solvekey::Symbol=:parametric,
                                    autodiff = :forward,
                                    algorithm=Optim.BFGS,
                                    algorithmkwargs=(), # add manifold to overwrite computed one
                                    options = Optim.Options(allow_f_increases=true,
                                                            time_limit = 100,
                                                            # show_trace = true,
                                                            # show_every = 1,
                                                            ))

  varIds = listVariables(fg)

  #TODO mabye remove sorting, just for convenience
  sort!(varIds, lt=natural_lt)
  separators = setdiff(varIds, frontals)

  varIds = [frontals; separators]

  flatvar = FlatVariables(fg, varIds)

  for vId in varIds
    p = getVariableSolverData(fg, vId, solvekey).val[1]
    flatvar[vId] =getCoordinates(getVariableType(fg,vId), p)
  end
  initValues = flatvar.X

  frontalsLength = sum(map(v->getDimension(getVariable(fg, v)), frontals))


  # build variables for frontals and seperators
  # fX = view(initValues, 1:frontalsLength)
  fX = initValues[1:frontalsLength]
  # sX = view(initValues, (frontalsLength+1):length(initValues))
  sX = initValues[frontalsLength+1:end]

  alg = algorithm(; algorithmkwargs...)
  # alg = algorithm(; algorithmkwargs...)
  cfd = calcFactorMahalanobisDict(fg)
  tdtotalCost = Optim.TwiceDifferentiable((x)->_totalCost(fg, cfd, flatvar, [x;sX]), fX, autodiff = autodiff)

  # result = Optim.optimize((x)->_totalCost(fg, flatvar, [x;sX]), fX, alg, options)
  result = Optim.optimize(tdtotalCost, fX, alg, options)

  rv = Optim.minimizer(result)
 
  H = Optim.hessian!(tdtotalCost, rv)

  Σ = pinv(H)

  d = OrderedDict{Symbol,NamedTuple{(:val, :cov),Tuple{Vector{Float64},Matrix{Float64}}}}()

  for key in frontals
    r = flatvar.idx[key]
    push!(d,key=>(val=rv[r],cov=Σ[r,r]))
  end

  return d, result, flatvar.idx, Σ
end

## ================================================================================================
## MixedCircular Manifold for Optim.jl
## ================================================================================================

"""
    MixedCircular
Mixed Circular Manifold. Simple manifold for circular and cartesian mixed for use in optim

DevNotes
- Consolidate around `ManifoldsBase.AbstractManifold` instead, with possible wrapper-type solution for `Optim.Manifold`
"""
struct MixedCircular <: Optim.Manifold
  isCircular::BitArray
end

function MixedCircular(fg::AbstractDFG, varIds::Vector{Symbol})
  circMask = Bool[]
  for k = varIds
    append!(circMask, convert(Tuple, getManifold(getVariableType(fg, k))) .== :Circular)
  end
  MixedCircular(circMask)
end

# https://github.com/JuliaNLSolvers/Optim.jl/blob/e439de4c997a727f3f724ae76da54b9cc08456b2/src/Manifolds.jl#L3
# retract!(m, x): map x back to a point on the manifold m
function Optim.retract!(c::MixedCircular, x)
  for (i,v) = enumerate(x)
    c.isCircular[i] && (x[i] = rem2pi(v, RoundNearest))
  end
  return x
end
# https://github.com/JuliaNLSolvers/Optim.jl/blob/e439de4c997a727f3f724ae76da54b9cc08456b2/src/Manifolds.jl#L2
# project_tangent!(m, g, x): project g on the tangent space to m at x
Optim.project_tangent!(S::MixedCircular,g,x) = g

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
  Aidx = IIF.collectIdx(flatvars,separators)
  Cidx = IIF.collectIdx(flatvars,frontals)

  #marginalize separators
  A = Σ[Aidx, Aidx]
  #marginalize frontals
  C = Σ[Cidx, Cidx]
  # cross
  B = Σ[Aidx, Cidx]


  Σₘ = deepcopy(A)
  if length(separators) == 0

    return (varlbl=Symbol[], μ=Float64[], Σ=Matrix{Float64}(undef,0,0))

  elseif length(separators) == 1

    # create messages
    return (varlbl = deepcopy(separators), μ = soldict[separators[1]].val, Σ = A)

  elseif length(separators) == 2
    A = Σₘ[1, 1]
    C = Σₘ[2, 2]
    B = Σₘ[1, 2]

    #calculate covariance between separators
    ΣA_B = A - B*inv(C)*B'
    # create messages
    m2lbl = deepcopy(separators)
    m2cov = isa(ΣA_B, Matrix) ? ΣA_B : fill(ΣA_B,1,1)
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
Add parametric solver to fg, batch solve using [`solveGraphParametric`](@ref) and update fg.
"""
function solveGraphParametric!(fg::AbstractDFG; init::Bool=true, kwargs...)
  if !(:parametric in fg.solverParams.algorithms)
    addParametricSolver!(fg; init=init)
  elseif init
    initParametricFrom!(fg)
  end

  vardict, result, varIds, Σ = solveGraphParametric(fg; kwargs...)

  updateParametricSolution!(fg, vardict)

  return vardict, result, varIds, Σ
end


"""
    $SIGNATURES
Initialize the parametric solver data from a different solution in `fromkey`.
"""
function initParametricFrom!(fg::AbstractDFG, fromkey::Symbol = :default; parkey::Symbol = :parametric, onepoint=false)
  if onepoint
    for v in getVariables(fg)
      fromvnd = getSolverData(v, fromkey)
      dims = getDimension(v)
      getSolverData(v, parkey).val[1] .= fromvnd.val[1]
      getSolverData(v, parkey).bw[1:dims,1:dims] .= LinearAlgebra.I(dims)
    end
  else
    for var in getVariables(fg)
        dims = getDimension(var)
        μ,Σ = calcMeanCovar(var, fromkey)
        getSolverData(var, parkey).val[1] .= μ
        getSolverData(var, parkey).bw[1:dims, 1:dims] .= Σ
    end
  end
end

"""
    $SIGNATURES
Add the parametric solveKey to all the variables in fg if it doesn't exists.
"""
function addParametricSolver!(fg; init=true)
  if !(:parametric in fg.solverParams.algorithms)
      push!(fg.solverParams.algorithms, :parametric)
      foreach(v->IIF.setDefaultNodeDataParametric!(v, getVariableType(v), initialized=false), getVariables(fg))
      if init
        initParametricFrom!(fg)
      end
  else
      error("parametric solvekey already exists")
  end
  nothing
end

#TODO Delete?
function updateVariablesFromParametricSolution!(fg::AbstractDFG, vardict)
  for (v,val) in vardict
    vnd = getVariableSolverData(fg, v, :parametric)
    vnd.val .= val.val
    if size(vnd.bw) != size(val.cov)
      vnd.bw = val.cov
    else
      vnd.bw .= val.cov
    end
  end
end

"""
    $SIGNATURES
Update the fg from solution in vardict and add MeanMaxPPE (all just mean). Usefull for plotting
"""
function updateParametricSolution!(sfg, vardict)

  for (v,val) in vardict
      vnd = getSolverData(getVariable(sfg, v), :parametric)
      # fill in the variable node data value
      p = getPoint(getVariableType(sfg, v), val.val)
      vnd.val[1] = p
      #calculate and fill in covariance
      vnd.bw = val.cov
      #fill in ppe as mean
      ppe = MeanMaxPPE(:parametric, val.val, val.val, val.val)
      getPPEDict(getVariable(sfg, v))[:parametric] = ppe
  end
end

function createMvNormal(val,cov)
    #TODO do something better for properly formed covariance, but for now just a hack...FIXME
    if all(diag(cov) .> 0.001) && isapprox(cov, transpose(cov), rtol=1e-4)
        return MvNormal(val,Symmetric(cov))
    else
        @error("Covariance matrix error", cov)
        # return nothing # FIXME, blanking nothing during #459 consolidation
        return MvNormal(val, ones(length(val)))
    end
end

function createMvNormal(v::DFGVariable, key=:parametric)
    if key == :parametric
        vnd = getSolverData(v, :parametric)
        dims = vnd.dims
        return createMvNormal(vnd.val[1:dims,1], vnd.bw[1:dims, 1:dims])
    else
        @warn "Trying MvNormal Fit, replace with PPE fits in future"
        return fit(MvNormal,getSolverData(v, key).val)
    end
end

## ================================================================================================
## Experimental specialized dispatch for Mixture
## ================================================================================================
# To sort out how to dispatch on specialized functions.
# related to #931 and #1069

struct MaxMixture
  p::Vector{Float64}
  # the chosen component to be used for the optimization
  choice::Base.RefValue{Int}
end

function getMeasurementParametric(s::Mixture{N,F,S,T}) where {N,F,S,T}
  meas = map(c->getMeasurementParametric(c)[1], values(s.components))
  iΣ = map(c->getMeasurementParametric(c)[2], values(s.components))
  return meas, iΣ
end

function _calcFactorMahalanobis(cfp, meas, iΣ, variables...)
  res = cfp.calcfactor!(meas, variables...)
  r = res' * iΣ * res
  return r
end

# DEV NOTE: function with other options including select once and use
# function (cfp::CalcFactorMahalanobis{<:CalcFactor, MaxMixture})(variables...)
#   if cfp.specialAlg.choice[] == 0
#     #calculate all mixture options
#     r = [_calcFactorMahalanobis(cfp, cfp.meas[i], cfp.iΣ[i], variables...) for i = 1:length(cfp.meas)]

#     p = cfp.specialAlg.p

#     k = size(cfp.iΣ[1], 2)
#     # α = 1 ./ sqrt.(2pi .* k .* det.(inv.(cfp.iΣ)))
#     α = sqrt.(det.(cfp.iΣ) ./ ((2pi)^k))

#     # mm, at = findmax(α .* p .* exp.(-0.5 .* r))
#     # mm = sum(α .* p .* exp.(-0.5 .* r) )
    
#     mm, at = findmin( 0.5 .* r .- log.(α .* p))
#     # mm = -log(sum(α .* p .* exp.(-0.5 .* r) ))
#     # return mm + maximum(log.(α .* p))
    
#     cfp.specialAlg.choice[] = at

#     return r[at] 

#   else
#     at = cfp.specialAlg.choice[]
#     return _calcFactorMahalanobis(cfp, cfp.meas[at], cfp.iΣ[at], variables...)
#   end

# end

# function (cfp::CalcFactorMahalanobis{<:CalcFactor, MaxMixture})(variables...)
  
#   r = [_calcFactorMahalanobis(cfp, cfp.meas[i], cfp.iΣ[i], variables...) for i = 1:length(cfp.meas)]

#   p = cfp.specialAlg.p

#   k = size(cfp.iΣ[1], 2)
#   # α = 1 ./ sqrt.(2pi .* k .* det.(inv.(cfp.iΣ)))
#   α = sqrt.(det.(cfp.iΣ) ./ ((2pi)^k))

#   mm, at = findmin(r .- log.(α .* p))
#   # mm = -log(sum(α .* p .* exp.(-0.5 .* r) ))
#   return mm + maximum(log.(α .* p))
    
# end


## ================================================================================================
## Experimental specialised dispatch for multihypo and nullhypo
## ================================================================================================
#TODO better dispatch

struct MaxMultihypo
  multihypo::Vector{Float64}
end
struct MaxNullhypo
  nullhypo::Float64
end

# function (cfp::CalcFactorMahalanobis{<:CalcFactor, MaxMultihypo})(X1, L1, L2)
#   mh = cfp.specialAlg.multihypo
#   @assert length(mh) == 3 "multihypo $mh  not supported with parametric, length should be 3"
#   @assert mh[1] == 0 "multihypo $mh  not supported with parametric, first should be 0"
  
#   #calculate both multihypo options
#   r1 = cfp(X1, L1)
#   r2 = cfp(X1, L2)
#   r = [r1, r2]

#   # hacky multihypo to start of with 
#   mm, at = findmin(r .* (1 .- mh[2:end]))
#   nat = at == 1 ? 1 : 2
#   k = length(X1)*one(r1) * 1e-3
#   return r[at] + r[nat]*k
  
# end

# function (cfp::CalcFactorMahalanobis{<:CalcFactor, MaxNullhypo})(X1, X2) 
#   nh = cfp.specialAlg.nullhypo
#   @assert nh > 0 "nullhypo $nh not as expected"
  
#   #calculate factor residual
#   res = cfp.calcfactor!(cfp.meas[1], X1, X2)
#   r1 =  res' * cfp.iΣ * res

#   # compare to uniform nullhypo
#   r2 = length(res)*one(r1)
#   r = [r1,r2]
#   mm, at = findmin(r .* [nh, (1-nh)])

#   residual = at == 1 ? r1 : r1*1e-3

#   return residual

#   # rand residual option
#   # idx = rand(Categorical([(1-nh), nh]))
#   # nh == 0.05 && cfp.varOrder==[:x1,:l1] && println("$idx -> $(r1.value), $r2")
#   # return r[idx] 

# end
