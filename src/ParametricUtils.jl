## ================================================================================================
## FlatVariables - used for packing variables for optimization
## ================================================================================================

struct FlatVariables{T<:Real}
  X::Vector{T}
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

function CalcFactorMahalanobis(fct::DFGFactor)
  cf = getFactorType(fct)
  varOrder = getVariableOrder(fct)
  
  _meas, _iΣ = getMeasurementParametric(cf)
  M = getManifold(getFactorType(fct))
  meas = typeof(_meas) <: Tuple ? _meas : (hat(M, Identity(M), _meas),)
  iΣ = typeof(_iΣ) <: Tuple ? _iΣ : (_iΣ,)

  # FIXME #1480, cache = preambleCache(dfg,vars,usrfnc)
  cache = nothing
  calcf = CalcFactor(getFactorMechanics(cf), _getFMdThread(fct), 0, 0, (), [], true, cache)
  
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

  return CalcFactorMahalanobis(calcf, varOrder, meas, iΣ, special)
end

# This is where the actual parametric calculation happens, CalcFactor equivalent for parametric
function (cfp::CalcFactorMahalanobis{CalcFactor{T, U, V, W, C}})(variables...) where {T<:AbstractFactor, U, V, W, C}
  # call the user function (be careful to call the new CalcFactor version only!!!)
  res = cfp.calcfactor!(cfp.meas[1], variables...)
  # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
  return res' * cfp.iΣ[1] * res
end

function (cfp::CalcFactorMahalanobis{CalcFactor{T, U, V, W, C}})(variables...) where {T<:AbstractPrior, U, V, W, C}
  # TODO should prior functions follow the same factor definition with a measurement on tangant?
  M = getManifold(cfp.calcfactor!.factor)
  #FIXME find a better fix for bitstypes
  ϵ = identity_element(M)
  ϵ = isbitstype(typeof(ϵ)) ? [ϵ] : ϵ
  m = exp(M, ϵ, cfp.meas[1])  

  res = cfp.calcfactor!(m, variables...)
  # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
  return res' * cfp.iΣ[1] * res
end

function (cfp::CalcFactorMahalanobis{CalcFactor{T, U, V, W, C}})(variables...) where {T<:Union{ManifoldFactor, ManifoldPrior} , U, V, W, C}
  # call the user function (be careful to call the new CalcFactor version only!!!)
  M = cfp.calcfactor!.factor.M
  X = cfp.calcfactor!(cfp.meas[1], variables...)
  # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
  # return mahalanobus_distance2(M, X, cfp.iΣ[1])
  
  #TODO do something about basis?
  # Xc = get_coordinates(M, variables[1], X, DefaultOrthogonalBasis())
  # Xc = get_coordinates(M, variables[1], X, DefaultOrthonormalBasis())
  Xc = vee(M, variables[1], X)
  return X, Xc' * cfp.iΣ[1] * Xc
end

function calcFactorMahalanobisDict(fg)
  calcFactors = Dict{Symbol, CalcFactorMahalanobis}()
  for fct in getFactors(fg)
    calcFactors[fct.label] = CalcFactorMahalanobis(fct)
  end
  return calcFactors
end

# new experimental Manifold ProductRepr version
function _totalCost(M, cfvec::Vector{<:CalcFactorMahalanobis}, labeldict, points)
  #
  obj = 0

  for cfp in cfvec

      varOrder = getindex.(Ref(labeldict), cfp.varOrder)

      factor_point_params = [points[M, i] for i in varOrder]

      # call the user function
      _,retval = cfp(factor_point_params...)

      # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
      obj += 1/2*retval
  end

  return obj
end

function _totalCost(fg, 
                    cfdict::Dict{Symbol, <:CalcFactorMahalanobis},
                    flatvar,
                    X )
  #
  obj = 0
  for (fid, cfp) in cfdict

    varOrder = cfp.varOrder

    Xparams = [getPoint(getVariableType(fg, varId), view(X, flatvar.idx[varId])) for varId in varOrder]

    # call the user function
    retval = cfp(Xparams...)

    # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
    obj += 1/2*retval
  end

  return obj
end


export solveGraphParametric

"""
    $SIGNATURES

Batch solve a Gaussian factor graph using Optim.jl. Parameters can be passed directly to optim.
Notes:
  - Only :Euclid and :Circular manifolds are currently supported, own manifold are supported with `algorithmkwargs` (code may need updating though)
"""
function solveGraphParametric(fg::AbstractDFG;
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

  #TODO mabye remove sorting, just for convenience
  sort!(varIds, lt=natural_lt)

  flatvar = FlatVariables(fg, varIds)

  for vId in varIds
    p = getVariableSolverData(fg, vId, solvekey).val[1]
    flatvar[vId] = getCoordinates(getVariableType(fg,vId), p)
  end

  initValues = flatvar.X


  alg = algorithm(; algorithmkwargs...)

  cfd = calcFactorMahalanobisDict(fg)
  tdtotalCost = Optim.TwiceDifferentiable((x)->_totalCost(fg, cfd, flatvar, x), initValues, autodiff = autodiff)

  result = Optim.optimize(tdtotalCost, initValues, alg, options)
  rv = Optim.minimizer(result)

  H = Optim.hessian!(tdtotalCost, rv)

  Σ = pinv(H)

  d = Dict{Symbol,NamedTuple{(:val, :cov),Tuple{Vector{Float64},Matrix{Float64}}}}()

  for key in varIds
    r = flatvar.idx[key]
    push!(d,key=>(val=rv[r],cov=Σ[r,r]))
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

  d = Dict{Symbol,NamedTuple{(:val, :cov),Tuple{Vector{Float64},Matrix{Float64}}}}()

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

function (cfp::CalcFactorMahalanobis{<:CalcFactor, MaxMixture})(variables...)
  
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

struct MaxMultihypo
  multihypo::Vector{Float64}
end
struct MaxNullhypo
  nullhypo::Float64
end

function (cfp::CalcFactorMahalanobis{<:CalcFactor, MaxMultihypo})(X1, L1, L2)
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

function (cfp::CalcFactorMahalanobis{<:CalcFactor, MaxNullhypo})(X1, X2) 
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
