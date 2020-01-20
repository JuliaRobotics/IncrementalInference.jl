
"""
    $SIGNATURES

Solve a Gaussian factor graph.
"""
function solveFactorGraphParametric(fg::AbstractDFG; solvekey::Symbol=:parametric)

  varIds = getVariableIds(fg)
  #TODO dimention di, its set to maximim and assumes all is the same
  vardims = map(v->getDimension(v), getVariables(fg))
  di = maximum(vardims)

  #TODO stoor dalk in kolomme
  initValues = zeros(di, length(varIds))
  varSymbols = Dict{Symbol, Int}()
  for (i,v) in enumerate(varIds)
    initValues[:,i] = getVal(getVariable(fg,v), 1; solveKey=solvekey)
    varSymbols[v] = i
  end

  #get all the cost functions of the factors.
  costfuns = []
  for fct in getFactors(fg)
    fac_cost_fx =  getFactorType(fct)

    varOrder = getVariableOrder(fct)
    # varsdata = [solverData(getVariable(fg, v), solvekey) for v in varOrder]
    #findfirst or dictionary lookup
    # idx = [findfirst(v->v==varId, varIds) for varId in varOrder]
    idx = [varSymbols[varId] for varId in varOrder]

    push!(costfuns, (fac_cost_fx, idx))

  end


  #build the cost function
  function totalCost(X)

    # dim = maximum(length.(X))
    dim = size(X)[1]

    res = zeros(dim)
    for (cf, idx) in costfuns
      Xparams = [X[:,i] for i in idx]
      res .+=  cf(Xparams...) #TODO splat performance?
    end
    # res[1]
    #TODO to get optim to work
    return norm(res)
  end

  result = optimize(totalCost, initValues, BFGS())
  rv = Optim.minimizer(result)

  d = Dict{Symbol,Vector{Float64}}()
  for (i,key) in enumerate(varIds)
    push!(d,key=>rv[:,i])
  end

  return d, result
end



#TODO maybe consolidate with solveFactorGraphParametric
function solveFrontalsParametric(fg::AbstractDFG, frontals::Vector{Symbol}; solvekey::Symbol=:parametric)

  varIds = getVariableIds(fg)
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
    # varsdata = [solverData(getVariable(fg, v), solvekey) for v in varOrder]
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
      nf = fit(Normal, solverData(var, fromkey).val)
      solverData(var, parkey).val[1,1] = nf.μ
      solverData(var, parkey).bw[1,1] = nf.σ
      # @show nf
      # m = var.estimateDict[:default].mean
  end
end
