
"""
    $SIGNATURES

Solve a Gaussian factor graph.
"""
function solveFactorGraphParametric!(fg::AbstractDFG; solvekey::Symbol=:parametric)

  varIds = getVariableIds(fg)
  #TODO dimention di
  di = 1
  #TODO stoor dalk in kolomme
  initValues = Array{Float64}(undef, di, length(varIds))
  varSymbols = Dict{Symbol, Int}()
  for (i,v) in enumerate(varIds)
    initValues[:,i] = getVal(getVariable(fg,v), 1; solveKey=solvekey)
    varSymbols[v] = i
  end
  # initValues = [getVal(getVariable(fg,v), 1; solveKey=solvekey)[1] for v in varIds]
  # TODO
  # do findfirst(v->v==:x0, varIds) lookup for or create dictionary
  # gebruik om trek te hou van indekse

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
  function totalCost(X)

    # dim = maximum(length.(X))
    dim = size(X)[1]

    res = zeros(dim)
    for (cf, idx) in costfuns
      res .+=  cf(X[:,idx]...) #TODO splat performance?
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
