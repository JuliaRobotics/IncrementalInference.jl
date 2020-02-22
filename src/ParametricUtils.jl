
"""
    $SIGNATURES

Solve a Gaussian factor graph.
"""
function solveFactorGraphParametric(fg::AbstractDFG; solvekey::Symbol=:parametric)

  varIds = listVariables(fg)
  #TODO dimention di, its set to maximim and assumes all is the same
  vardims = map(v->getDimension(v), getVariables(fg))

  shapes = [ArrayShape{Real}(d) for d in vardims]
  varshapes = NamedTupleShape(NamedTuple{Tuple(varIds)}(shapes))

  initValues = zeros(totalndof(varshapes))

  function totalCost(X)

    shapedX = varshapes(X)
    res::Float64 = 0
    for fct in getFactors(fg)

      cf = getFactorType(fct)
      varOrder = getVariableOrder(fct)

      Xparams = [collect(getproperty(shapedX, varId)) for varId in varOrder]

      retval = cf(Xparams...)
      # @assert retVal |> typeof == Float64
      res += 1/2*retval # 1/2*log(1/(  sqrt(det(Σ)*(2pi)^k) ))  ## k = dim(μ)
    end
    return res

  end

  tdtotalCost = TwiceDifferentiable(totalCost, initValues)
  result = optimize(tdtotalCost, initValues, BFGS())
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
function solveFrontalsParametric(fg::AbstractDFG, frontals::Vector{Symbol}; solvekey::Symbol=:parametric)

  varIds = listVariables(fg)
  separators = setdiff(varIds, frontals)

  #TODO dimention di, its set to maximim and assumes all is the same
  vardims = map(v->getDimension(v), getVariables(fg))
  di = maximum(vardims)

  initValues = Array{Float64}(undef, di, length(varIds))
  varSymbols = Dict{Symbol, Int}()
  i = 1

  # pack frontals first and then seperators in the same array
  for v in frontals
    initValues[:,i] = getVal(getVariable(fg,v), 1; solveKey=solvekey)
    varSymbols[v] = i
    i +=1
  end

  for v in separators
    initValues[:,i] = getVal(getVariable(fg,v), 1; solveKey=solvekey)
    varSymbols[v] = i
    i +=1
  end

  #get all the cost functions of the factors.
  costfuns = []
  for fct in getFactors(fg)
    fac_cost_fx =  getFactorType(fct)

    varOrder = getVariableOrder(fct)
    # varsdata = [getSolverData(getVariable(fg, v), solvekey) for v in varOrder]
    #TODO findfirst or dictionary lookup
    # idx = [findfirst(v->v==varId, varIds) for varId in varOrder]
    idx = [varSymbols[varId] for varId in varOrder]

    push!(costfuns, (fac_cost_fx, idx))

  end


  #build the cost function
  function totalCost(F, S)
    X = [F S]
    # dim = maximum(length.(X))
    dim = size(X)[1]

    res = zeros(dim)
    for (cf, idx) in costfuns
      Xparams = [X[:,i] for i in idx]
      res .+=  cf(Xparams...) #TODO splat performance?
      # res .+=  cf(X[:,idx]...) #TODO splat performance?
    end
    # res[1]
    #TODO to get optim to work
    return norm(res)
  end

  # build variables for frontals and seperators
  fX = initValues[:,1:length(frontals)]
  sX = initValues[:,length(frontals)+1:end]

  result = optimize(x->totalCost(x,sX), fX, BFGS())
  rv = Optim.minimizer(result)


  d = Dict{Symbol,Vector{Float64}}()
  for (i,key) in enumerate(frontals)
    push!(d,key=>rv[:,i])
  end

  return d, result

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
