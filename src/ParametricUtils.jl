
"""
    $SIGNATURES

Solve a Gaussian factor graph.
"""
function solveFactorGraphParametric(fg::AbstractDFG;
                                    solvekey::Symbol=:parametric,
                                    autodiff = :finite,
                                    algorithm=BFGS(),
                                    options = Optim.Options())

  varIds = listVariables(fg)

  #TODO remove sorting, just for convenience
  sort!(varIds, lt=natural_lt)


  vardims = map(v->getDimension(getVariable(fg, v)), varIds)

  shapes = [ArrayShape{Real}(d) for d in vardims]
  varshapes = NamedTupleShape(NamedTuple{Tuple(varIds)}(shapes))

  initValues = zeros(totalndof(varshapes))
  # initValues = Vector{Float64}(undef, totalndof(varshapes))

  shapedinitValues = varshapes(initValues)
  #populate initial values from current fg values
  for vId in varIds
    getproperty(shapedinitValues,vId) .= getVariableSolverData(fg, vId, solvekey).val[:,1]
  end

  function totalCost(X)

    shapedX = varshapes(X)
    res = 0
    for fct in getFactors(fg)

      cf = getFactorType(fct)
      varOrder = getVariableOrder(fct)

      Xparams = [getproperty(shapedX, varId) for varId in varOrder]

      retval = cf(Xparams...)
      # @assert retVal |> typeof == Float64
      res += 1/2*retval # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
    end
    return res

  end

  tdtotalCost = TwiceDifferentiable(totalCost, initValues, autodiff = autodiff)
  result = optimize(tdtotalCost, initValues, algorithm, options)
  rv = Optim.minimizer(result)

  H = hessian!(tdtotalCost, rv)

  Σ = pinv(H)

  d = Dict{Symbol,NamedTuple{(:val, :cov),Tuple{Vector{Float64},Matrix{Float64}}}}()
  rvshaped = varshapes(rv)

  for key in varIds
    s = getproperty(varshapes, key)
    r = range(s.offset+1, length=s.len)
    push!(d,key=>(val=getproperty(rvshaped,key),cov=Σ[r,r]))
  end

  return d, result, varshapes, Σ
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
                                    autodiff = :finite,
                                    algorithm=BFGS(),
                                    options = Optim.Options())

  varIds = listVariables(fg)
  separators = setdiff(varIds, frontals)


  varIdsFS = [frontals; separators]

  vardims = map(v->getDimension(getVariable(fg, v)), varIdsFS)

  #TODO look at ConstValueShape
  shapes = [ArrayShape{Real}(d) for d in vardims]

  frontalsLength = sum(map(v->getDimension(getVariable(fg, v)), frontals))


  varshapes = NamedTupleShape(NamedTuple{Tuple(varIdsFS)}(shapes))

  initValues = zeros(totalndof(varshapes))
  # initValues = Array{Float64}(undef, di, length(varIds))
  # initValues = Vector{Float64}(undef, totalndof(varshapes))

  shapedinitValues = varshapes(initValues)
  #populate initial values from current fg values
  for vId in varIdsFS
    getproperty(shapedinitValues,vId) .= getVariableSolverData(fg, vId, solvekey).val[:,1]
  end



  #build the cost function
  function totalCost(X)

    shapedX = varshapes(X)
    res = 0
    for fct in getFactors(fg)

      cf = getFactorType(fct)
      varOrder = getVariableOrder(fct)

      Xparams = [getproperty(shapedX, varId) for varId in varOrder]

      retval = cf(Xparams...)
      # @assert retVal |> typeof == Float64
      res += 1/2*retval # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
    end
    return res

  end

  # build variables for frontals and seperators
  fX = initValues[1:frontalsLength]
  sX = initValues[frontalsLength+1:end]

  tdtotalCost = TwiceDifferentiable(totalCost, initValues, autodiff = autodiff)
  # result = optimize(x->totalCost(x,sX), fX, BFGS())

  totalCost([fX;sX])

  result = optimize(x->totalCost([x;sX]), fX, algorithm, options)

  rv = Optim.minimizer(result)

  H = hessian!(tdtotalCost, [rv; sX])

  Σ = pinv(H)

  d = Dict{Symbol,NamedTuple{(:val, :cov),Tuple{Vector{Float64},Matrix{Float64}}}}()
  rvshaped = varshapes([rv; sX])

  for key in varIds
    s = getproperty(varshapes, key)
    r = range(s.offset+1, length=s.len)
    push!(d,key=>(val=getproperty(rvshaped,key),cov=Σ[r,r]))
  end

  return d, result, varshapes, Σ
end


"""
    $SIGNATURES
Get the indexes for labels in shape varShape
"""
function collectIdx(varShape, labels)
  idx = Int[]
  for lbl in labels
    append!(idx, collect(getproperty(varShape,lbl).offset.+(1:getproperty(varShape,lbl).len)))
  end
  return idx
end

"""
    $SIGNATURES
Get the indexes for labels in shape varShape
"""
function calculateCoBeliefMessage(soldict, Σ, varIdsShape, separators, frontals)
  Aidx = IIF.collectIdx(varIdsShape,separators)
  Cidx = IIF.collectIdx(varIdsShape,frontals)

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
function initParametricFrom(fg::AbstractDFG, fromkey::Symbol = :default; parkey::Symbol = :parametric)
  for var in getVariables(fg)
      #TODO only supports Normal now
      # expand to MvNormal
      nf = fit(Normal, getSolverData(var, fromkey).val)
      getSolverData(var, parkey).val[1,1] = nf.μ
      getSolverData(var, parkey).bw[1,1] = nf.σ
      # @show nf
      # m = var.estimateDict[:default].mean
  end
end
