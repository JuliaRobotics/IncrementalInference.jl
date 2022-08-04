# ================================================================================================
## FactorOperationalMemory for parametric, TODO move back to FactorOperationalMemory.jl
## ================================================================================================

abstract type AbstractMaxMixtureSolver end

"""
$TYPEDEF

Internal parametric extension to [`CalcFactor`](@ref) used for buffering measurement and calculating Mahalanobis distance

Related

[`CalcFactor`](@ref)
"""
struct CalcFactorMahalanobis{N, S<:Union{Nothing,AbstractMaxMixtureSolver}}
  calcfactor!::CalcFactor
  varOrder::Vector{Symbol}
  meas::NTuple{N, <:AbstractArray}
  iΣ::NTuple{N, Matrix{Float64}}
  specialAlg::S
end
# struct CalcFactorMahalanobis{CF<:CalcFactor, S<:Union{Nothing,AbstractMaxMixtureSolver}, N}
#   calcfactor!::CF
#   varOrder::Vector{Symbol}
#   meas::NTuple{N, <:AbstractArray}
#   iΣ::NTuple{N, Matrix{Float64}}
#   specialAlg::S
# end


# ================================================================================================
## FlatVariables - used for packing variables for optimization
## ================================================================================================

struct FlatVariables{T<:Real}
  X::Vector{T}
  idx::OrderedDict{Symbol, UnitRange{Int}}
end

function FlatVariables(fg::AbstractDFG, varIds::Vector{Symbol})

  index = 1
  idx = OrderedDict{Symbol, UnitRange{Int}}()
  for vid = varIds
    v = getVariable(fg, vid)
    dims = getDimension(v)
    idx[vid] = index:(index+dims-1)
    index += dims
  end
  return FlatVariables(Vector{Float64}(undef, index-1), idx)
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
  else
    @warn "getMeasurementParametric falls back to using field `.Z` by default. Extend it for more complex factors." 
    error("getMeasurementParametric(::$(typeof(s))) not defined, please add it, or use non-parametric, or open an issue for help.")
  end
  
  return getMeasurementParametric(Z)
end



## ================================================================================================
## Parametric solve with Mahalanobis distance - CalcFactor
## ================================================================================================

#TODO maybe remove with Mixture rework see #1504
getFactorMechanics(f::AbstractFactor) = f
getFactorMechanics(f::Mixture) = f.mechanics

function CalcFactorMahalanobis(fg, fct::DFGFactor)
  fac_func = getFactorType(fct)
  varOrder = getVariableOrder(fct)
  
  _meas, _iΣ = getMeasurementParametric(fac_func)
  M = getManifold(getFactorType(fct))

  ϵ = getPointIdentity(M)
  
  if typeof(_meas) <: Tuple
    _measX = map(m->hat(M, ϵ, m), _meas)
  # TODO perhaps better consolidate manifold prior 
  elseif fac_func isa ManifoldPrior
    _measX = (_meas,)
  else
    _measX = (hat(M, ϵ, _meas),)
  end

  meas = fac_func isa AbstractPrior ? map(X->exp(M, ϵ, X), _measX) : _measX

  iΣ = typeof(_iΣ) <: Tuple ? _iΣ : (_iΣ,)

  cache = preambleCache(fg, getVariable.(fg, varOrder), getFactorType(fct))

  # calcf = CalcFactor(getFactorMechanics(fac_func), nothing, 0, 0, nothing, nothing, true, nothing)
  calcf = CalcFactor(getFactorMechanics(fac_func), nothing, 0, 0, nothing, nothing, true, cache)
  #FactorMetadata is not type stable
  # calcf = CalcFactor(getFactorMechanics(fac_func), _getFMdThread(fct), 0, 0, nothing, nothing, true, cache)
  
  multihypo = getSolverData(fct).multihypo
  nullhypo = getSolverData(fct).nullhypo

  # FIXME, type instability, use dispatch instead of if-else
  if length(multihypo) > 0
    special = MaxMultihypo(multihypo)
  elseif nullhypo > 0 
    special = MaxNullhypo(nullhypo)
  elseif fac_func isa Mixture
    special = MaxMixture(fac_func.diversity.p, Ref(0))
  else
    special = nothing
  end

  return CalcFactorMahalanobis(calcf, varOrder, meas, iΣ, special)
end


# This is where the actual parametric calculation happens, CalcFactor equivalent for parametric
function (cfp::CalcFactorMahalanobis{1,Nothing})(variables...) # AbstractArray{T} where T <: Real
  # call the user function 
  res = cfp.calcfactor!(cfp.meas[1], variables...)
  # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
  return res' * cfp.iΣ[1] * res
end


function calcFactorMahalanobisDict(fg)
  calcFactors = OrderedDict{Symbol, CalcFactorMahalanobis}()
  for fct in getFactors(fg)
    calcFactors[fct.label] = CalcFactorMahalanobis(fg, fct)
  end
  return calcFactors
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
## GraphSolveStructures
## ================================================================================================


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
  #TODO tuple or vector?
  # vartypes = tuple(keys(typedict)...)
  vartypes = collect(keys(typedict))
  return vartypes, typedict, alltypes
end

function buildGraphSolveManifold(fg::AbstractDFG)
  vartypes, vartypecount, vartypeslist = getVariableTypesCount(fg)

  PMs = map(vartypes) do vartype
      N = vartypecount[vartype] 
      G = getManifold(vartype)
      NPowerManifold(G, N)
      # PowerManifold(G, NestedReplacingPowerRepresentation(), N)
      # PowerManifold(G, NestedPowerRepresentation(), N) #TODO investigate as it does not converge
  end
  M = ProductManifold(PMs...)    
  return M, vartypes, vartypeslist
end

struct GraphSolveBuffers{T<:Real, U}
  ϵ::U
  p::U
  X::U
  Xc::Vector{T}
end

function GraphSolveBuffers(@nospecialize(M), ::Type{T}) where T
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
  varOrderDict::OrderedDict{Symbol, Tuple{Int, Vararg{Int,}}}
end


function GraphSolveContainer(fg)

  M, varTypes, varTypesIds = buildGraphSolveManifold(fg)
  varTypesIndexes = ArrayPartition(values(varTypesIds)...)
  buffs = OrderedDict{DataType,GraphSolveBuffers}()
  cfd = calcFactorMahalanobisDict(fg)

  varOrderDict = OrderedDict{Symbol, Tuple{Int, Vararg{Int,}}}()
  for (fid, cfp) in cfd 
    varOrder = cfp.varOrder
    var_idx = map(varOrder) do v
      findfirst(==(v), varTypesIndexes)
    end
    varOrderDict[fid] = tuple(var_idx...)
  end

  return GraphSolveContainer(M, buffs, varTypes, varTypesIds, cfd, varOrderDict)
end

function getGraphSolveCache!(gsc::GraphSolveContainer, ::Type{T}) where T<:Real
  cache = gsc.buffers
  M = gsc.M
  val = get!(cache, T) do 
    @debug "cache miss, cacheing" T
    GraphSolveBuffers(M, T)
  end
  return val
end

function _toPoints2!(@nospecialize(M::AbstractManifold), buffs::GraphSolveBuffers{T,U}, Xc::Vector{T}) where {T,U}
  ϵ = buffs.ϵ
  p = buffs.p
  X = buffs.X
  get_vector!(M, X, ϵ, Xc, DefaultOrthogonalBasis())
  exp!(M, p, ϵ, X)
  return p::U
end

cost_cfp(cfp::CalcFactorMahalanobis, p, vi::NTuple{1,Int}) = cfp(p[vi[1]])
cost_cfp(cfp::CalcFactorMahalanobis, p, vi::NTuple{2,Int}) = cfp(p[vi[1]], p[vi[2]])
cost_cfp(cfp::CalcFactorMahalanobis, p, vi::NTuple{3,Int}) = cfp(p[vi[1]], p[vi[2]], p[vi[3]])

# the cost function
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

#NOTE this only works with a product of power manifolds
function getComponentsCovar(@nospecialize(PM::ProductManifold), Σ::AbstractMatrix)

  dims = manifold_dimension.(PM.manifolds)
  dim_ranges = Manifolds._get_dim_ranges(dims)

  subsigmas = map(zip(dim_ranges, PM.manifolds)) do v
    r = v[1]
    M = v[2]
    _getComponentsCovar(M, view(Σ, r, r)) 
  end

  return ArrayPartition(subsigmas...)

end

function _getComponentsCovar(@nospecialize(PM::PowerManifold), Σ::AbstractMatrix)
  M = PM.manifold
  dim = manifold_dimension(M)
  subsigmas = map(Manifolds.get_iterator(PM)) do i
    r = (i-1)*dim+1:i*dim
    Σ[r, r]
  end

  return subsigmas
end

function solveGraphParametric(fg::AbstractDFG;
                              computeCovariance::Bool = true,
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
  initValues .+= randn(length(Xc))*0.0001

  #optim setup and solve
  alg = algorithm(; algorithmkwargs...)
  tdtotalCost = Optim.TwiceDifferentiable(gsc, initValues, autodiff = autodiff)

  result = Optim.optimize(tdtotalCost, initValues, alg, options)
  rv = Optim.minimizer(result)

  # optionally compute hessian for covariance
  Σ = if computeCovariance
    H = Optim.hessian!(tdtotalCost, rv)
    pinv(H)
  else
    N = length(initValues)
    zeros(N,N)
  end

  #TODO better return 

  #get point (p) values form results
  get_vector!(M, X, ϵ, rv, DefaultOrthogonalBasis())
  exp!(M, p, ϵ, X)

  #extract covariances from result
  # sigmas = getComponentsCovar(M, Σ)

  # d = OrderedDict{Symbol,NamedTuple{(:val, :cov),Tuple{Vector{Float64},Matrix{Float64}}}}()
  d = OrderedDict{Symbol,NamedTuple{(:val, :cov),Tuple{AbstractArray,Matrix{Float64}}}}()

  varIds = vcat(values(gsc.varTypesIds)...)
  varIdDict = FlatVariables(fg, varIds).idx
  for (i,key) in enumerate(varIds)
    r = varIdDict[key]
    push!(d,key=>(val=p[i], cov=Σ[r,r]))
    # push!(d,key=>(val=p[i], cov=sigmas[i]))
  end


  return (opti=d, stat=result, varIds=varIdDict, Σ=Σ)
end


## Original
# ==============================

function _totalCost(fg, 
                    cfdict::OrderedDict{Symbol, <:CalcFactorMahalanobis},
                    flatvar,
                    Xc )
  #
  obj = zero(eltype(Xc))
  for (fid, cfp) in cfdict

    varOrder = cfp.varOrder

    Xparams = [getPoint(getVariableType(fg, varId), view(Xc, flatvar.idx[varId])) for varId in varOrder]

    # call the user function
    retval = cfp(Xparams...)

    # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
    obj += 1/2*retval
  end

  return obj
end


#TODO maybe consolidate with solveGraphParametric
#TODO WIP
```
    $SIGNATURES
Solve for frontal values only with values in seprarators fixed
```
function solveConditionalsParametric(fg::AbstractDFG,
                                    frontals::Vector{Symbol},
                                    separators::Vector{Symbol} = setdiff(listVariables(fg), frontals);
                                    solvekey::Symbol=:parametric,
                                    autodiff = :forward,
                                    algorithm=Optim.BFGS,
                                    algorithmkwargs=(), # add manifold to overwrite computed one
                                    options = Optim.Options(allow_f_increases=true,
                                                            time_limit = 100,
                                                            # show_trace = true,
                                                            # show_every = 1,
                                                            ))

  varIds = [frontals; separators]

  sfg = issetequal(varIds, listVariables(fg)) ? fg : buildSubgraph(fg, varIds, 1)

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
  cfd = calcFactorMahalanobisDict(sfg)
  tdtotalCost = Optim.TwiceDifferentiable((x)->_totalCost(fg, cfd, flatvar, [x;sX]), fX, autodiff = autodiff)

  # result = Optim.optimize((x)->_totalCost(fg, flatvar, [x;sX]), fX, alg, options)
  result = Optim.optimize(tdtotalCost, fX, alg, options)

  rv = Optim.minimizer(result)
 
  H = Optim.hessian!(tdtotalCost, rv)

  Σ = pinv(H)

  d = OrderedDict{Symbol,NamedTuple{(:val, :cov),Tuple{AbstractArray,Matrix{Float64}}}}()

  for key in frontals
    r = flatvar.idx[key]
    p = getPoint(getVariableType(fg, key), rv[r])
    push!(d,key=>(val=p,cov=Σ[r,r]))
  end
  
  return (opti=d, stat=result, varIds=flatvar.idx, Σ=Σ)
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

"""
    $SIGNATURES
Update the fg from solution in vardict and add MeanMaxPPE (all just mean). Usefull for plotting
"""
function updateParametricSolution!(sfg, vardict)

  for (v,val) in vardict
      vnd = getSolverData(getVariable(sfg, v), :parametric)
      # fill in the variable node data value
      vnd.val[1] = val.val
      #calculate and fill in covariance
      vnd.bw = val.cov
      #fill in ppe as mean
      Xc = getCoordinates(getVariableType(sfg, v), val.val)

      ppe = MeanMaxPPE(:parametric, Xc, Xc, Xc)
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


function autoinitParametric!(dfg::AbstractDFG,
                             initme::Symbol;
                             reinit::Bool=false)
  #
  solveKey = :parametric

  xi = getVariable(dfg, initme)
  vnd = getSolverData(xi, solveKey)
  # don't initialize a variable more than once
  if reinit || !isInitialized(xi, solveKey)

    # frontals - initme
    # separators - inifrom

    initfrom = ls2(dfg, initme)
    filter!(initfrom) do vl
      isInitialized(dfg, vl, solveKey)
    end

    vardict, result, flatvars, Σ = solveConditionalsParametric(dfg, [initme], initfrom)

    val,cov = vardict[initme]

    # fill in the variable node data value
    vnd.val[1] = val
    #calculate and fill in covariance
    vnd.bw = cov

    vnd.initialized = true
    #fill in ppe as mean
    Xc = getCoordinates(getVariableType(xi), val)
    ppe = MeanMaxPPE(:parametric, Xc, Xc, Xc)
    getPPEDict(xi)[:parametric] = ppe

    # updateVariableSolverData!(dfg, xi, solveKey, true; warn_if_absent=false)    
    # updateVariableSolverData!(dfg, xi.label, getSolverData(xi, solveKey), :graphinit, true, Symbol[]; warn_if_absent=false)
  else
    result = nothing
  end

  return result#isInitialized(xi, solveKey)
end

## ================================================================================================
## Experimental specialized dispatch for Mixture
## ================================================================================================
# To sort out how to dispatch on specialized functions.
# related to #931 and #1069

struct MaxMixture <: AbstractMaxMixtureSolver
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

function (cfp::CalcFactorMahalanobis{N, MaxMixture})(variables...) where N
  
  r = [_calcFactorMahalanobis(cfp, cfp.meas[i], cfp.iΣ[i], variables...) for i = 1:length(cfp.meas)]

  p = cfp.specialAlg.p

  k = size(cfp.iΣ[1], 2)
  # α = 1 ./ sqrt.(2pi .* k .* det.(inv.(cfp.iΣ)))
  α = sqrt.(det.(cfp.iΣ) ./ ((2pi)^k))

  mm, at = findmin(r .- log.(α .* p))
  # mm = -log(sum(α .* p .* exp.(-0.5 .* r) ))
  return mm + maximum(log.(α .* p))
    
end


## ================================================================================================
## Experimental specialised dispatch for multihypo and nullhypo
## ================================================================================================
#TODO better dispatch

struct MaxMultihypo <: AbstractMaxMixtureSolver
  multihypo::Vector{Float64}
end
struct MaxNullhypo <: AbstractMaxMixtureSolver
  nullhypo::Float64
end

function (cfp::CalcFactorMahalanobis{N, MaxMultihypo})(X1, L1, L2) where N
  mh = cfp.specialAlg.multihypo
  @assert length(mh) == 3 "multihypo $mh  not supported with parametric, length should be 3"
  @assert mh[1] == 0 "multihypo $mh  not supported with parametric, first should be 0"
  
  #calculate both multihypo options
  r1 = cfp(X1, L1)
  r2 = cfp(X1, L2)
  r = [r1, r2]

  # hacky multihypo to start of with 
  mm, at = findmin(r .* (1 .- mh[2:end]))
  nat = at == 1 ? 1 : 2
  k = length(X1)*one(r1) * 1e-3
  return r[at] + r[nat]*k
  
end

function (cfp::CalcFactorMahalanobis{N, MaxNullhypo})(X1, X2) where N
  nh = cfp.specialAlg.nullhypo
  @assert nh > 0 "nullhypo $nh not as expected"
  
  #calculate factor residual
  res = cfp.calcfactor!(cfp.meas[1], X1, X2)
  r1 =  res' * cfp.iΣ * res

  # compare to uniform nullhypo
  r2 = length(res)*one(r1)
  r = [r1,r2]
  mm, at = findmin(r .* [nh, (1-nh)])

  residual = at == 1 ? r1 : r1*1e-3

  return residual

  # rand residual option
  # idx = rand(Categorical([(1-nh), nh]))
  # nh == 0.05 && cfp.varOrder==[:x1,:l1] && println("$idx -> $(r1.value), $r2")
  # return r[idx] 

end
