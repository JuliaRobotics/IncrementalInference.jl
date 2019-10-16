
function compareField(Allc, Bllc, syms)::Bool
  return eval(:($Allc.$syms == $Bllc.$syms))
end

"""
    $(SIGNATURES)

Compare the all fields of T that are not in `skip` for objects `Al` and `Bl`.

TODO > add to func_ref.md
"""
function compareFields(Al::T,
                       Bl::T;
                       show::Bool=true,
                       skip::Vector{Symbol}=Symbol[]  )::Bool where {T}
  TP = true
  fields = fieldnames(T)
  for field in fields
    if (field in skip)
      continue
    end
    tp = compareField(Al, Bl, field)
    show ? println("$tp : $field") : nothing
    TP = TP && tp
  end
  return TP
end

function compareFields(Al::T,
                       Bl::T;
                       show::Bool=true,
                       skip::Vector{Symbol}=Symbol[]  )::Bool where {T <: Union{Number, AbstractString}}
  #
  return Al == Bl
end

function compareAll(Al::T,
                    Bl::T;
                    show::Bool=true,
                    skip::Vector{Symbol}=Symbol[]  )::Bool where {T <: Union{AbstractString,Symbol}}
  #
  return Al == Bl
end

function compareAll(Al::T,
                    Bl::T;
                    show::Bool=true,
                    skip::Vector{Symbol}=Symbol[]  )::Bool where {T <: Union{Array{<:Number}, Number}}
  #
  if length(Al) != length(Bl)
    return false
  end
  return norm(Al - Bl) < 1e-6
end

function compareAll(Al::T,
                    Bl::T;
                    show::Bool=true,
                    skip::Vector{Symbol}=Symbol[]  )::Bool where {T <: Array}
  #
  if length(Al) != length(Bl)
    return false
  end
  TP = true
  for i in 1:length(Al)
    TP = TP && compareAll(Al[i],Bl[i], show=false)
  end
  return TP
end


"""
    $(SIGNATURES)

Recursively compare the all fields of T that are not in `skip` for objects `Al` and `Bl`.

TODO > add to func_ref.md
"""
function compareAll(Al::T,
                    Bl::T;
                    show::Bool=true,
                    skip::Vector{Symbol}=Symbol[]  )::Bool where {T <: Tuple}
  #
  TP = true
  TP = TP && length(Al) == length(Bl)
  for i in 1:length(Al)
    TP &= compareAll(Al[i], Bl[i], show=show, skip=skip)
  end
  return TP
end

function compareAll(Al::T,
                    Bl::T;
                    show::Bool=true,
                    skip::Vector{Symbol}=Symbol[]  )::Bool where {T <: Dict}
  #
  TP = true
  TP = TP && length(Al) == length(Bl)
  !TP ? (return false) : nothing
  for (id, val) in Al
    if Symbol(id) in skip
      continue
    end
    TP = TP && compareAll(val, Bl[id], show=show, skip=skip)
  end
  return TP
end

# TODO: KEEP
function compareAll(Al::T1, Bl::T2; show::Bool=true, skip::Vector{Symbol}=Symbol[])::Bool where {T1 <: Union{SingleThreaded, MultiThreaded}, T2 <: Union{SingleThreaded, MultiThreaded}}
  return T1 == T2
end

function compareAll(Al::T, Bl::T; show::Bool=true, skip::Vector{Symbol}=Symbol[])::Bool where T
  @show T
  @show Al
  @show Bl
  TP = compareFields(Al, Bl, show=show, skip=skip)
  !TP ? (return false;) : nothing
  for field in fieldnames(T)
    println("field: $field")
    if field in skip
      continue
    end
    Ad = eval(:($Al.$field))
    Bd = eval(:($Bl.$field))
    TP = TP && compareAll(Ad, Bd, show=show, skip=skip)
    println(string(TP))
  end
  return TP
end

# function compareAll(Al::T,
#                     Bl::T;
#                     show::Bool=true,
#                     skip::Vector{Symbol}=Symbol[]  )::Bool where {T <: Dict}
#   #
#   TP = true
#   TP = TP && length(Al) == length(Bl)
#   !TP ? (return false) : nothing
#   for (id, val) in Al
#     if Symbol(id) in skip
#       continue
#     end
#     compareAll(val, Bl[id], show=show, skip=skip)
#   end
#   return true
# end

# TODO: KEEP
function compare(p1::BallTreeDensity, p2::BallTreeDensity)::Bool
  return compareAll(p1.bt,p2.bt, skip=[:calcStatsHandle; :data]) &&
         compareAll(p1,p2, skip=[:calcStatsHandle; :bt])
end

"""
    $SIGNATURES

Compare that all fields are the same in a `::FactorGraph` variable.
"""
function compareVariable(A::DFGVariable,
                         B::DFGVariable;
                         skip::Vector{Symbol}=Symbol[],
                         show::Bool=true,
                         skipsamples::Bool=true  )::Bool
  #
  skiplist = union([:attributes;:solverDataDict;:_internalId],skip)
  TP = compareAll(A, B, skip=skiplist, show=show)
  varskiplist = skipsamples ? [:val; :bw] : Symbol[]
  skiplist = union([:softtype;],varskiplist)
  union!(skiplist, skip)
  TP = TP && compareAll(A.solverDataDict, B.solverDataDict, skip=skiplist, show=show)

  Ad = solverData(A)
  Bd = solverData(B)

  # TP = TP && compareAll(A.attributes, B.attributes, skip=[:softtype;], show=show)
  varskiplist = union(varskiplist, [:softtype;:_internalId])
  union!(varskiplist, skip)
  TP = TP && compareAll(Ad, Bd, skip=varskiplist, show=show)
  TP = TP && typeof(Ad.softtype) == typeof(Bd.softtype)
  TP = TP && compareAll(Ad.softtype, Bd.softtype, show=show, skip=skip)
  return TP
end

function compareAllSpecial(A::T1,
                           B::T2;
                           skip=Symbol[],
                           show::Bool=true) where {T1 <: GenericFunctionNodeData, T2 <: GenericFunctionNodeData}
  if T1 != T2
    return false
  else
    return compareAll(A, B, skip=skip, show=show)
  end
end

# TODO: KEEP
function compareAllSpecial(A::T1, B::T2;
                    skip=Symbol[], show::Bool=true) where {T1 <: CommonConvWrapper, T2 <: CommonConvWrapper}
  #
  if T1 != T2
    return false
  else
    return compareAll(A, B, skip=skip, show=show)
  end
end

"""
    $SIGNATURES

Compare that all fields are the same in a `::FactorGraph` factor.
"""
function compareFactor(A::DFGFactor,
                       B::DFGFactor;
                       show::Bool=true,
                       skip::Vector{Symbol}=Symbol[],
                       skipsamples::Bool=true,
                       skipcompute::Bool=true  )
  #
  @info show
  TP =  compareAll(A, B, skip=union([:attributes;:data;:_variableOrderSymbols;:_internalId],skip), show=show)
  # TP = TP & compareAll(A.attributes, B.attributes, skip=[:data;], show=show)
  TP = TP & compareAllSpecial(solverData(A), solverData(B), skip=union([:fnc;:_internalId], skip), show=show)
  TP = TP & compareAllSpecial(solverData(A).fnc, solverData(B).fnc, skip=union([:cpt;:measurement;:params;:varidx;:threadmodel], skip), show=show)
  TP = TP & (skipsamples || compareAll(solverData(A).fnc.measurement, solverData(B).fnc.measurement, show=show, skip=skip))
  TP = TP & (skipcompute || compareAll(solverData(A).fnc.params, solverData(B).fnc.params, show=show, skip=skip))
  TP = TP & (skipcompute || compareAll(solverData(A).fnc.varidx, solverData(B).fnc.varidx, show=show, skip=skip))

  return TP
end



"""
    $SIGNATURES

Compare all variables in both `::FactorGraph`s A and B.

Notes
- A and B should all the same variables and factors.

Related:

`compareFactorGraphs`, `compareSimilarVariables`, `compareVariable`, `ls`
"""
function compareAllVariables(fgA::G1,
                             fgB::G2;
                             skip::Vector{Symbol}=Symbol[],
                             show::Bool=true,
                             skipsamples::Bool=true )::Bool where {G1 <: AbstractDFG, G2 <: AbstractDFG}
  # get all the variables in A or B
  xlA =  getVariableIds(fgA)
  xlB =  getVariableIds(fgB)
  vars = union(xlA, xlB)

  # compare all variables exist in both A and B
  TP = length(xlA) == length(xlB)
  for xla in xlA
    TP &= xla in xlB
  end
  # slightly redundant, but repeating opposite direction anyway
  for xlb in xlB
    TP &= xlb in xlA
  end

  # compare each variable is the same in both A and B
  for var in vars
    TP = TP && compareVariable(DFG.getVariable(fgA, var), DFG.getVariable(fgB, var), skipsamples=skipsamples, skip=skip)
  end

  # return comparison result
  return TP
end

"""
    $SIGNATURES

Compare similar labels between `::FactorGraph`s A and B.

Notes
- At least one variable label should exist in both A and B.

Related:

`compareFactorGraphs`, `compareAllVariables`, `compareSimilarFactors`, `compareVariable`, `ls`.
"""
function compareSimilarVariables(fgA::G1,
                                 fgB::G2;
                                 skip::Vector{Symbol}=Symbol[],
                                 show::Bool=true,
                                 skipsamples::Bool=true )::Bool where {G1 <: AbstractDFG, G2 <: AbstractDFG}
  #
  xlA = getVariableIds(fgA)
  xlB = getVariableIds(fgB)

  # find common variables
  xlAB = intersect(xlA, xlB)
  TP = length(xlAB) > 0

  # compare the common set
  for var in xlAB
    @info var
    TP &= compareVariable(DFG.getVariable(fgA, var), DFG.getVariable(fgB, var), skipsamples=skipsamples, skip=skip)
  end

  # return comparison result
  return TP
end

"""
    $SIGNATURES

Determine if and compare `fgS::FactorGraph` is a subset with similar content to `fgA`.

Notes
- `fgS` ⊆ `fgA`.

Related:

`compareFactorGraphs`, `compareSimilarVariables`, `compareSimilarFactors`, `ls`.
"""
function compareSubsetFactorGraph(fgS::FactorGraph, fgA::FactorGraph)
  error("not implemented yet")
  return false
end

"""
    $SIGNATURES

Compare similar factors between `::FactorGraph`s A and B.

Related:

`compareFactorGraphs`, `compareSimilarVariables`, `compareAllVariables`, `ls`.
"""
function compareSimilarFactors(fgA::G1,
                               fgB::G2;
                               skipsamples::Bool=true,
                               skipcompute::Bool=true,
                               show::Bool=true  )::Bool where {G1 <: AbstractDFG, G2 <: AbstractDFG}
  #
  xlA = getFactorIds(fgA)
  xlB = getFactorIds(fgB)

  # find common variables
  xlAB = intersect(xlA, xlB)
  TP = length(xlAB) > 0

  # compare the common set
  for var in xlAB
    TP = TP && compareFactor(DFG.getFactor(fgA, var), getFactor(fgB, var), skipsamples=skipsamples, skipcompute=skipcompute, show=show)
  end

  # return comparison result
  return TP
end

"""
    $SIGNATURES

Compare and return if two factor graph objects are the same, by comparing similar variables and factors.

Notes:
- Default items to skip with `skipsamples`, `skipcompute`.
- User defined fields to skip can be specified with `skip::Vector{Symbol}`.

Related:

`compareSimilarVariables`, `compareSimilarFactors`, `compareAllVariables`, `ls`.
"""
function compareFactorGraphs(fgA::G1,
                             fgB::G2;
                             skipsamples::Bool=true,
                             skipcompute::Bool=true,
                             skip::Vector{Symbol}=Symbol[],
                             show::Bool=true  )::Bool where {G1 <: AbstractDFG, G2 <: AbstractDFG}
  #
  skiplist = Symbol[:g;:bn;:IDs;:fIDs;:id;:nodeIDs;:factorIDs;:fifo;:solverParams]
  skiplist = union(skiplist, skip)
  @warn "compareFactorGraphs will skip comparisons on: $skiplist"

  TP = compareAll(fgA, fgB, skip=skiplist, show=show)
  TP = TP && compareSimilarVariables(fgA, fgB, skipsamples=skipsamples, show=show, skip=skiplist )
  TP = TP && compareSimilarFactors(fgA, fgB, skipsamples=skipsamples, skipcompute=skipcompute, show=show )
  TP = TP && compareAll(fgA.solverParams, fgB.solverParams, skip=skiplist)

  return TP
end
