
export
  getCliqVariableInferredPercent,
  getCliqVariableMoreInitDims,
  getVariablePossibleDim


## Description
# Dim              -- Native dimension of Variable/Factor
# Possible         -- Maximum dimension possible in variable
# Solvable         -- number of dimensions which can be resolved from current state
# Inferred         -- Current inferered dimension available
# suffix-Fraction  -- Report as percentage fraction
# XPercentage      -- Report ratio of X over Possible


## Major objectives
#
# getCliqVariableInferredPercent
# getCliqVariableMoreInitDims


## Variables

"""
    $SIGNATURES

Return the number of dimensions this variable vertex `var` contains.

Related

getVariableInferredDim, getVariableInferredDimFraction
"""
getVariableDim(vard::VariableNodeData)::Int = getSofttype(vard).dims
getVariableDim(var::DFGVariable)::Int = getVariableDim(solverData(var))

"""
    $SIGNATURES

Return the number of projected dimensions into a variable during inference.

Notes
- `saturate` clamps return value to no greater than variable dimension

Related

getVariableDim, getVariableInferredDimFraction, getVariableInferredDim, getVariableDim
"""
getVariableInferredDim(vard::VariableNodeData, saturate::Bool=false) = saturate && getVariableDim(vard) < vard.inferdim ? getVariableDim(vard) : vard.inferdim
getVariableInferredDim(var::DFGVariable, saturate::Bool=false) = getVariableInferredDim(solverData(var), saturate)
function getVariableInferredDim(fg::G, varid::Symbol, saturate::Bool=false) where G <: AbstractDFG
  getVariableInferredDim(getVariable(fg, varid), saturate)
end


getVariableInferredDimFraction(vard::VariableNodeData, saturate::Bool=false)::Float64 = getVariableInferredDim(vard, saturate) / getVariableDim(vard)
getVariableInferredDimFraction(var::DFGVariable, saturate::Bool=false)::Float64 = getVariableInferredDim(solverData(var), saturate)
function getVariableInferredDimFraction(dfg::G, varid::Symbol, saturate::Bool=false)::Float64 where G <: AbstractDFG
  getVariableInferredDimFraction(getVariable(dfg, varid), saturate)
end


## Factors

"""
    $SIGNATURES

Return the number of dimensions this factor vertex `fc` influences.
"""
getFactorDim(fcd::GenericFunctionNodeData)::Float64 = isa(fcd.fnc.usrfnc!, MsgPrior) ? fcd.fnc.usrfnc!.inferdim : Float64(fcd.fnc.zDim)
getFactorDim(fc::DFGFactor)::Int = getFactorDim(solverData(fc))
function getFactorDim(fg::G, fctid::Symbol)::Int where G <: AbstractDFG
  getFactorDim(getFactor(fg, fctid))
end

"""
   $SIGNATURES

Return the sum of factor dimensions connected to variable as per the factor graph `fg`.

Related

getFactorSolvableDim, getVariableDim, getVariableInferredDim, getFactorDim, isCliqFullDim
"""
function getVariablePossibleDim(fg::G, var::DFGVariable, fcts::Vector{Symbol}=ls(fg, var.label))::Float64 where G <: AbstractDFG
  alldims = 0.0
  for fc in fcts
    alldims += getFactorDim(fg, fc)
  end
  return alldims
end
function getVariablePossibleDim(fg::G,
                                varid::Symbol,
                                fcts::Vector{Symbol}=ls(fg, varid)  )::Float64 where G <: AbstractDFG
  #
  getVariablePossibleDim(fg, getVariable(fg, varid))
end

"""
    $SIGNATURES

Return the total inferred dimension available for variable from factor based on current inferred status of other connected variables.

Notes
- Accumulate the factor dimension fractions:  Sum [0..1]*zdim
- Variable dimenion fractions are inferdim / vardim
- Variable dimension are saturated at vardim for the calculating solve dimensions

Related

getVariablePossibleDim, getVariableDim, getVariableInferredDim, getFactorDim, getFactorSolvableDim, isCliqFullDim
"""
function getFactorInferFraction(dfg::G,
                                idfct::Symbol,
                                varid::Symbol,
                                saturate::Bool=false  )::Float64 where G <: AbstractDFG
  # get all other variables
  allvars = lsf(dfg, idfct)
  lievars = setdiff(allvars, [varid;])

  # get all other var dimensions with saturation
  len = length(lievars)
  fracs = map(lv->getVariableInferredDimFraction(dfg, lv, true), lievars)

  if length(fracs) == 0
    return 0.0
  end

  # the dimension of leave one out variables dictate if this factor can prodive full information on leave out variable.
  return cumprod(fracs)[end]
end


"""
    $SIGNATURES

Return the total inferred/solvable dimension available for variable based on current inferred status of other factor connected variables.

Notes
- Accumulate the factor dimension fractions:  Sum [0..1]*zdim
- Variable dimenion fractions are inferdim / vardim
- Variable dimension are saturated at vardim for the calculating solve dimensions

Related

getVariablePossibleDim, getVariableDim, getVariableInferredDim, getFactorDim, getFactorInferFraction, isCliqFullDim
"""
function getFactorSolvableDimFraction(dfg::G,
                                      idfct::Symbol,
                                      varid::Symbol,
                                      saturate::Bool=false  )::Float64 where G <: AbstractDFG
  #
  # get all other variables
  allvars = lsf(dfg, idfct)

  # prior/unary
  if length(allvars) == 1
    return 1.0
  end

  # general case
  lievars = setdiff(allvars, [varid;])

  # get all other var dimensions with saturation
  len = length(lievars)
  fracs = map(lv->getVariableInferredDimFraction(dfg, lv, true), lievars)

  # the dimension of leave one out variables dictate if this factor can prodive full information on leave out   variable.
  return cumprod(fracs)[end]
end
function getFactorSolvableDimFraction(dfg::G,
                                      fct::DFGFactor,
                                      varid::Symbol,
                                      saturate::Bool=false  )::Float64 where G <: AbstractDFG
  #
  getFactorSolvableDimFraction(dfg,fct.label,varid,saturate)
end

function getFactorSolvableDim(dfg::G,
                              idfct::Symbol,
                              varid::Symbol,
                              saturate::Bool=false  )::Float64 where G <: AbstractDFG
  #
  return getFactorSolvableDimFraction(dfg,idfct,varid,saturate)*getFactorDim(dfg, idfct)
end
function getFactorSolvableDim(dfg::G,
                              fct::DFGFactor,
                              varid::Symbol,
                              saturate::Bool=false  )::Float64 where G <: AbstractDFG
  #
  return getFactorSolvableDimFraction(dfg,fct,varid,saturate)*getFactorDim(fct)
end

"""
    $SIGNATURES

Return the total solvable dimension for each variable in the factor graph `dfg`.

Notes
- "Project" the solved dimension from other variables through connected factors onto each variable.
"""
function getVariableSolvableDim(dfg::G, varid::Symbol, fcts::Vector{Symbol}=ls(dfg, varid)) where G <: AbstractDFG

  sd = 0.0
  for fc in fcts
    sd += getFactorSolvableDim(dfg,fc,varid)
  end
  return sd
end


## Combined Variables and Factors

"""
    $SIGNATURES

Return the current dimensionality of solve for each variable in a clique.
"""
function getCliqVariableInferDims(dfg::G,
                                  cliq::Graphs.ExVertex,
                                  saturate::Bool=true,
                                  fraction::Bool=true  )::Dict{Symbol,Float64} where G <: AbstractDFG
  #
  # which variables
  varids = getCliqAllVarIds(cliq)

  # and what inferred dimension in this dfg
  retd = Dict{Symbol,Float64}()
  for varid in varids
    retd[varid] = getVariableInferredDim(dfg, varid)
  end

  return retd
end

# """
#     $SIGNATURES
#
# Return the directly achievable dimensionality of solve for each variable in a clique.
#
# Related
#
# getFactorSolvableDim
# """
# function getCliqVarPossibleDim(dfg::G,
#                                     cliq::Graphs.ExVertex,
#                                     saturate::Bool=true,
#                                     fraction::Bool=true  )::Dict{Symbol, Float64}
#   #
#   # variables and factors associated with this clique
#   vars = getCliqAllVarIds(cliq)
#   # fcts = getCliqAllFactIds(cliq)
#   # rows = length(fcts)
#   cols = length(vars)
#
#   # for output result
#   dict = Dict{Symbol,Float64}()
#
#   for j in 1:cols
#     dict[vars[j]] += getVariablePossibleDim(dfg, vars[j])
#   end
#
# end


"""
    $SIGNATURES

Return dictionary of clique variables and percentage of inference completion for each.

Notes
- Completion means (relative to clique subgraph) ratio of inferred dimension over possible solve dimension.

Related

getCliqVariableMoreInitDims
"""
function getCliqVariableInferredPercent(dfg::G, cliq::Graphs.ExVertex) where G <: AbstractDFG

  # cliq variables factors
  vars = getCliqAllVarIds(cliq)
  # fcts = getCliqAllFactIds(cliq)
  # rows = length(fcts)
  # nvars = length(vars)

  # for output result
  dict = Dict{Symbol,Float64}()

  # possible variable infer dim
  # current variable infer dim
  # calculate ratios in [0,1]
  for var in getCliqAllVarIds(cliq)# 1:nvars
      dict[var] = getVariableInferredDim(dfg, var)
      dict[var] /= getVariablePossibleDim(dfg, var)
  end

  # return dict with result
  return dict
end


"""
    $SIGNATURES

Return a dictionary with the number of immediately additionally available inference
dimensions on each variable in a clique.

Related

getCliqVariableInferredPercent
"""
function getCliqVariableMoreInitDims(dfg::G,
                                     cliq::Graphs.ExVertex  ) where G <: AbstractDFG
  #
  # cliq variables factors
  vars = getCliqAllVarIds(cliq)

  # for output result
  dict = Dict{Symbol,Float64}()

  # possible variable infer dim
  # current variable infer dim
  for vari in vars
    dict[vari] = getVariableSolvableDim(dfg, vari)
    dict[vari] -= getVariableInferredDim(dfg, vari)
  end

  # return dict with result
  return dict
end


"""
    $SIGNATURES

Return true if the variables solve dimension is equal to the sum of connected factor dimensions.

Related

getVariableInferredDimFraction, getVariableDim, getVariableInferredDim, getVariablePossibleDim
"""
function isCliqFullDim(fg::G,
                       cliq::Graphs.ExVertex )::Bool where G <: AbstractDFG
  #
  # get various variable percentages
  red = getCliqVariableInferredPercent(fg, cliq)

  # if all variables are solved to their full potential
  return abs(sum(collect(values(red))) - length(red)) < 1e10
end



"""
    $SIGNATURES

Return a vector of clique index `::Int` for descending order solvable dimensions of each clique.

Notes
- Uses the number of inferable/solvable dimension information in descending order.
- Uses function fetchCliqSolvableDims/getCliqVariableMoreInitDims to retrieved cached values of new solvable/inferable dimensions.

Related

fetchCliqSolvableDims, getCliqVariableMoreInitDims, getSubFgPriorityInitOrder
"""
function getCliqSiblingsPriorityInitOrder(tree::BayesTree,
                                          prnt::Graphs.ExVertex,
                                          logger=ConsoleLogger() )::Vector{Int}
  #
  sibs = getChildren(tree, prnt)
  len = length(sibs)
  tdims = Vector{Int}(undef, len)
  sidx = Vector{Int}(undef, len)
  for idx in 1:len
    cliqd = getData(sibs[idx])
    with_logger(logger) do
      @info "getCliqSiblingsPriorityInitOrder, idx=$idx of $len, $(cliqd.frontalIDs[1]) length solvableDims=$(length(cliqd.solvableDims.data))"
    end
    flush(logger.stream)
    sidims = fetchCliqSolvableDims(sibs[idx])
    sidx[idx] = sibs[idx].index
    tdims[idx] = sum(collect(values(sidims)))
    with_logger(logger) do
      @info "getCliqSiblingsPriorityInitOrder, finished idx=$idx of $len, length solvableDims=$(length(cliqd.solvableDims.data))"
    end
    flush(logger.stream)
  end
  p = sortperm(tdims, rev=true)
  with_logger(logger) do
    @info "getCliqSiblingsPriorityInitOrder, done p=$p"
  end
  return sidx[p]
end


"""
    $SIGNATURES

Return a vector of variable labels `::Symbol` for descending order solvable dimensions of each clique.

Notes
- EXPERIMENTAL, NOT EXPORTED
- Uses the number of inferable/solvable dimension information in descending order.
- Uses function getVariableSolvableDim to retrieved cached values of new solvable/inferable dimensions.

Related

getVariableSolvableDim, getCliqSiblingsPriorityInitOrder
"""
function getSubFgPriorityInitOrder(sfg::G, logger=ConsoleLogger()) where G <: AbstractDFG

    vars = ls(sfg)
    len = length(vars)
    tdims = Vector{Int}(undef, len)
    for idx in 1:len
      tdims[idx] = getVariableSolvableDim(sfg, vars[idx])
    end
    p = sortperm(tdims, rev=true)
    with_logger(logger) do
      @info "getSubFgPriorityInitOrder -- ordered vars=$(vars[p])"
      @info "getSubFgPriorityInitOrder -- ordered tdims=$(tdims[p])"
    end
    return vars[p]
end
