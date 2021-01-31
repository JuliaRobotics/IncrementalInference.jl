

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

"""
    $SIGNATURES

Returns the parametric measurement for a factor as a tuple (measurement, inverse covariance) for parametric inference (assumign Gaussian).
Defaults to find the parametric measurement at field `Z` followed by `z`.  
  
Notes
- Users should overload this method should their factor not default to `.Z<:ParametricType` or `.z<:ParametricType`
"""
function getParametricMeasurement end

function getParametricMeasurement(Z)
  error("$Z not supported, please use non-parametric or open an issue if it should be")
end

function getParametricMeasurement(Z::Normal)
  meas = mean(Z)
  iσ = 1/std(Z)^2
  return [meas], reshape([iσ],1,1)
end

function getParametricMeasurement(Z::MvNormal)
  meas = mean(Z)
  iΣ = invcov(Z)
  return meas, iΣ
end

function getParametricMeasurement(s::FunctorInferenceType)
  if hasfield(typeof(s), :Z)
    Z = s.Z
  elseif hasfield(typeof(s), :z)
    Z = s.z
  else
    error("$s not supported, please use non-parametric or open an issue if it should be")
  end
  return getParametricMeasurement(Z)
end


"""
    $TYPEDEF

Internal parametric extension to [`CalcFactor`](@ref) 
"""
struct CalcFactorMahalanobis{CF}
  calcfactor!::CF
  varOrder::Vector{Symbol}
  meas::Vector{Float64}
  iΣ::Matrix{Float64}
end

function CalcFactorMahalanobis(fct::DFGFactor)
  cf = getFactorType(fct)
  varOrder = getVariableOrder(fct)
  meas, iΣ = getParametricMeasurement(cf)
  calcf = CalcFactor(cf, _getFMdThread(fct), 0, 0, (), [])
  return CalcFactorMahalanobis(calcf, varOrder, meas, iΣ)
end

# pass in residual for consolidation with nonparametric
# userdata is now at `cfp.cf.cachedata`
function (cfp::CalcFactorMahalanobis)(variables...)
  # call the user function (be careful to call the new CalcFactor version only!!!)
  res = Vector{eltype(variables[1])}(undef, length(cfp.meas))
  iΣ = cfp.iΣ
  cfp.calcfactor!(res, cfp.meas, variables...)
  
  # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
  return 0.5 * (res' * iΣ * res)
end

function calcFactorMahalanobisDict(fg)
  calcFactors = Dict{Symbol, CalcFactorMahalanobis}()
  for fct in getFactors(fg)
    calcFactors[fct.label] = CalcFactorMahalanobis(fct)
  end
  return calcFactors
end

function _totalCost(cfdict::Dict{Symbol, CalcFactorMahalanobis}, 
                    flatvar, 
                    X )
  #
  obj = 0
  for (fid, cfp) in cfdict

    varOrder = cfp.varOrder
    
    Xparams = [view(X, flatvar.idx[varId]) for varId in varOrder]

    # call the user function
    retval = cfp(Xparams...) 

    # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
    obj += 1/2*retval 
  end

  return obj
end


# build the cost function
function _totalCost(fg::AbstractDFG, 
                    flatvar, 
                    X )
  #
  obj = 0
  for fct in getFactors(fg)

    cf = getFactorType(fct)
    varOrder = getVariableOrder(fct)
    
    Xparams = [view(X, flatvar.idx[varId]) for varId in varOrder]

    retval = cf(Xparams...)

    # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
    obj += 1/2*retval 
  end

  return obj
end


"""
    $SIGNATURES

Solve a Gaussian factor graph.
"""
function solveFactorGraphParametric(fg::AbstractDFG;
                                    useCalcFactor::Bool=false, #TODO dev param will be removed
                                    solvekey::Symbol=:parametric,
                                    autodiff = :forward,
                                    algorithm=BFGS,
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
  varIds = listVariables(fg)

  #TODO mabye remove sorting, just for convenience
  sort!(varIds, lt=natural_lt)

  flatvar = FlatVariables(fg, varIds)

  for vId in varIds
    flatvar[vId] = getVariableSolverData(fg, vId, solvekey).val[:,1]
  end

  initValues = flatvar.X


  mc_mani = IIF.MixedCircular(fg, varIds)
  alg = algorithm(;manifold=mc_mani, algorithmkwargs...)
  # alg = algorithm(; algorithmkwargs...)

  if useCalcFactor
    cfd = calcFactorMahalanobisDict(fg)
    tdtotalCost = TwiceDifferentiable((x)->_totalCost(cfd, flatvar, x), initValues, autodiff = autodiff)
  else
    tdtotalCost = TwiceDifferentiable((x)->_totalCost(fg, flatvar, x), initValues, autodiff = autodiff)
  end

  result = optimize(tdtotalCost, initValues, alg, options)
  rv = Optim.minimizer(result)

  H = hessian!(tdtotalCost, rv)

  Σ = pinv(H)

  d = Dict{Symbol,NamedTuple{(:val, :cov),Tuple{Vector{Float64},Matrix{Float64}}}}()

  for key in varIds
    r = flatvar.idx[key]
    push!(d,key=>(val=rv[r],cov=Σ[r,r]))
  end

  return d, result, flatvar.idx, Σ
end

#TODO maybe consolidate with solveFactorGraphParametric
#TODO WIP
```
    $SIGNATURES
Solve for frontal values only with values in seprarators fixed
```
function solveConditionalsParametric(fg::AbstractDFG,
                                    frontals::Vector{Symbol};
                                    solvekey::Symbol=:parametric,
                                    autodiff = :forward,
                                    algorithm=BFGS,
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
    flatvar[vId] = getVariableSolverData(fg, vId, solvekey).val[:,1]
  end
  initValues = flatvar.X

  frontalsLength = sum(map(v->getDimension(getVariable(fg, v)), frontals))


  # build variables for frontals and seperators
  # fX = view(initValues, 1:frontalsLength)
  fX = initValues[1:frontalsLength]
  # sX = view(initValues, (frontalsLength+1):length(initValues))
  sX = initValues[frontalsLength+1:end]

  mc_mani = MixedCircular(fg, varIds)
  alg = algorithm(;manifold=mc_mani, algorithmkwargs...)

  tdtotalCost = TwiceDifferentiable((x)->_totalCost(fg, flatvar,x), initValues, autodiff = autodiff)

  result = optimize((x)->_totalCost(fg, flatvar, [x;sX]), fX, alg, options)
  # result = optimize(x->totalCost([x;sX]), fX, alg, options)

  rv = Optim.minimizer(result)

  H = hessian!(tdtotalCost, [rv; sX])

  Σ = pinv(H)

  d = Dict{Symbol,NamedTuple{(:val, :cov),Tuple{Vector{Float64},Matrix{Float64}}}}()

  for key in frontals
    r = flatvar.idx[key]
    push!(d,key=>(val=rv[r],cov=Σ[r,r]))
  end

  return d, result, flatvar.idx, Σ
end

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

"""
    $SIGNATURES

Initialize the parametric solver data from a different solution in `fromkey`.
"""
function initParametricFrom!(fg::AbstractDFG, fromkey::Symbol = :default; parkey::Symbol = :parametric)
  for var in getVariables(fg)
      #TODO only supports Normal now
      # expand to MvNormal
      fromvnd = getSolverData(var, fromkey)
      if fromvnd.dims == 1
        nf = fit(Normal, fromvnd.val)
        getSolverData(var, parkey).val[1,1] = nf.μ
        getSolverData(var, parkey).bw[1,1] = nf.σ
        # @show nf
        # m = var.estimateDict[:default].mean
      else
        #FIXME circular will not work correctly with MvNormal
        nf = fit(MvNormal, fromvnd.val)
        getSolverData(var, parkey).val[1:fromvnd.dims] .= mean(nf)[1:fromvnd.dims]
        getSolverData(var, parkey).bw = cov(nf)
      end
  end
end

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
    MixedCircular
Mixed Circular Manifold. Simple manifold for circular and cartesian mixed for use in optim
"""
struct MixedCircular <: Optim.Manifold
  isCircular::BitArray
end

function MixedCircular(fg::AbstractDFG, varIds::Vector{Symbol})
  circMask = Bool[]
  for k = varIds
    append!(circMask, getVariableType(fg, k) |> getManifolds .== :Circular)
  end
  MixedCircular(circMask)
end

function Optim.retract!(c::MixedCircular, x)
  for (i,v) = enumerate(x)
    c.isCircular[i] && (x[i] = rem2pi(v, RoundNearest))
  end
  return x
end
Optim.project_tangent!(S::MixedCircular,g,x) = g


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




## SANDBOX of usefull development functions to be cleaned up
"""
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
Update the fg from solution in vardict and add MeanMaxPPE (all just mean). Usefull for plotting
"""
function updateParametricSolution(sfg, vardict)

  for (v,val) in vardict
      vnd = getSolverData(getVariable(sfg, v), :parametric)
      # fill in the variable node data value
      vnd.val .= val.val
      #calculate and fill in covariance
      vnd.bw = val.cov
      #fill in ppe as mean
      ppe = MeanMaxPPE(:parametric, val.val, val.val, val.val) 
      getPPEDict(getVariable(sfg, v))[:parametric] = ppe
  end
end