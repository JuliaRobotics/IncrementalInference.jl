# Factor Graph OS type utilities
#  IIF methods should direclty detect extended types from user import
# of convert in their namespace



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

  Ad = getData(A)
  Bd = getData(B)

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
                       skipsamples::Bool=true,
                       skipcompute::Bool=true  )
  #
  @info show
  TP =  compareAll(A, B, skip=[:attributes;:data;:_variableOrderSymbols;:_internalId], show=show)
  # TP = TP & compareAll(A.attributes, B.attributes, skip=[:data;], show=show)
  TP = TP & compareAllSpecial(getData(A), getData(B), skip=[:fnc;:_internalId], show=show)
  TP = TP & compareAllSpecial(getData(A).fnc, getData(B).fnc, skip=[:cpt;:measurement;:params;:varidx;:threadmodel], show=show)
  TP = TP & (skipsamples || compareAll(getData(A).fnc.measurement, getData(B).fnc.measurement, show=show))
  TP = TP & (skipcompute || compareAll(getData(A).fnc.params, getData(B).fnc.params, show=show))
  TP = TP & (skipcompute || compareAll(getData(A).fnc.varidx, getData(B).fnc.varidx, show=show))

  return TP
end
  # Ad = getData(A)
  # Bd = getData(B)
  # TP =  compareAll(A, B, skip=[:attributes;:data], show=show)
  # TP &= compareAll(A.attributes, B.attributes, skip=[:data;], show=show)
  # TP &= compareAllSpecial(getData(A).fnc, getData(B).fnc, skip=[:cpt;], show=show)
  # TP &= compareAll(getData(A).fnc.cpt, getData(B).fnc.cpt, show=show)


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
function compareSubsetFactorGraph(fgS::FactorGraph, fgA::FactorGraph; api::DataLayerAPI=localapi)
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
  skiplist = Symbol[:g;:bn;:IDs;:fIDs;:id;:nodeIDs;:factorIDs;:fifo; :solverParams]
  skiplist = union(skiplist, skip)
  @warn "compareFactorGraphs will skip comparisons on: $skiplist"

  TP = compareAll(fgA, fgB, skip=skiplist, show=show)
  TP = TP && compareSimilarVariables(fgA, fgB, skipsamples=skipsamples, show=show, skip=skiplist )
  TP = TP && compareSimilarFactors(fgA, fgB, skipsamples=skipsamples, skipcompute=skipcompute, show=show )
  TP = TP && compareAll(fgA.solverParams, fgB.solverParams, skip=skiplist)

  return TP
end

"""
    $(SIGNATURES)

Save mostly complete Factor Graph type by converting complicated FunctionNodeData
types to 'Packed' types using user supplied converters. Ground truth can also be
saved and recovered by the associated loadjld(file="tempfg.jld2") method.

Notes:
- Must use `.jld2` since Julia 1.0 (previous version was deprecated).
"""
function savejld(fgl::G;
                 file::AbstractString="tempfg.jld2"  ) where G <: AbstractDFG
  # fgs = encodefg(fgl)
  @save file fgl
  return file
end



"""
    $(SIGNATURES)

Opposite of savejld(fg, gt=gt, file="tempfg.jl") to load data from file. This function
uses the unpacking converters for converting all PackedInferenceType to FunctorInferenceType.
"""
function loadjld(;file::AbstractString="tempfg.jld2")
  fgd = @load file fgl
  return fgd
end

"""
    $(SIGNATURES)

Test if all elements of the string is a number:  Ex, "123" is true, "1_2" is false.
"""
allnums(str::S) where {S <: AbstractString} = occursin(Regex(string(["[0-9]" for j in 1:length(str)]...)), str)
# occursin(r"_+|,+|-+", node_idx)

isnestednum(str::S; delim='_') where {S <: AbstractString} = occursin(Regex("[0-9]+$(delim)[0-9]+"), str)

function sortnestedperm(strs::Vector{<:AbstractString}; delim='_')
  str12 = split.(strs, delim)
  sp1 = sortperm(parse.(Int,getindex.(str12,2)))
  sp2 = sortperm(parse.(Int,getindex.(str12,1)[sp1]))
  return sp1[sp2]
end


"""
    $(SIGNATURES)

Return array of all variable nodes connected to the last `n` many poses (`:x*`).

Example:

```julia
# Shallow copy the tail end of poses from a factor graph `fg1`
vars = lsRear(fg1, 5)
fg1_r5 = subgraphFromVerts(fg1, vars)
```
"""
function lsRear(fgl::FactorGraph, n::Int=1)
  lasts = ls(fgl)[1][(end-n):end]
  syms = ls(fgl, lasts)
  union(lsf.(fgl, syms)[:]...)
end

"""
    $SIGNATURES

Return `::Bool` on whether `fg::FactorGraph` has orphaned nodes or graph fragments.
"""
hasOrphans(fg::FactorGraph) = sum(length.(ls.(fg, [ls(fg)[1];ls(fg)[2]])) .== 0) > 0

"""
    $SIGNATURES

Return `Vector{Symbol}` of landmarks attached to vertex vsym in `fgl::FactorGraph`.
"""
function landmarks(fgl::FactorGraph, vsym::Symbol)
  fsyms = ls(fgl, vsym)
  lms = Symbol[]
  for fs in fsyms
    for varv = lsf(fgl, fs)
      if string(varv)[1] == 'l'
        push!(lms, varv)
      end
    end
  end
  lms
end



function evalLikelihood(fg::FactorGraph, sym::Symbol, point::Vector{Float64})
  p = getVertKDE(fg, sym)
  Ndim(p) == length(point) ? nothing : error("point (dim=$(length(point))) must have same dimension as belief (dim=$(Ndim(p)))")
  evaluateDualTree(p, reshape(point,:,1))[1]
end

# Evaluate the likelihood of an Array{2} of points on the marginal belief of some variable
# note the dimensions must match
function evalLikelihood(fg::FactorGraph, sym::Symbol, points::Array{Float64,2})
  p = getVertKDE(fg, sym)
  Ndim(p) == size(points,1) ? nothing : error("points (dim=$(size(points,1))) must have same dimension as belief (dim=$(Ndim(p)))")
  evaluateDualTree(p, (points))
end



function setThreadModel!(fgl::FactorGraph;model=IncrementalInference.SingleThreaded)
  for (key, id) in fgl.fIDs
    getData(getVert(fgl, key,nt=:fnc)).fnc.threadmodel = model
  end
  nothing
end

function _evalType(pt::String)::Type
    try
        getfield(Main, Symbol(pt))
    catch ex
        io = IOBuffer()
        showerror(io, ex, catch_backtrace())
        err = String(take!(io))
        error("_evalType: Unable to locate factor/distribution type '$pt' in main context (e.g. do a using on all your factor libraries). Please check that this factor type is loaded into main. Stack trace = $err")
    end
end

"""
    $(SIGNATURES)

Print the maximum point values form all variables approximate marginals in the factor graph.
The full marginal can be recovered for example `X0 = getVertKDE(fg, :x0)`.
"""
function printgraphmax(fgl::FactorGraph)
    verts = union(ls(fgl)...)
    map(v -> println("$v : $(getKDEMax(getVertKDE(fgl, v)))"), verts);
end


"""
    $SIGNATURES

Return interger index of desired variable element.

Example
-------
```julia
pp = RoME.Point2()
getIdx(pp, :posY) # = 2
```

Internal Notes
--------------
- uses number i < 100 for index number, and
- uses +100 offsets to track the minibatch number of the requested dimension
"""
function getIdx(pp::T, sym::Symbol, i::Int=0)::Tuple{Int, Int} where {T <: Tuple}
  # i > 99 ? error("stop") : nothing
  i-=100
  for p in pp
    i,j = getIdx(p, sym, i)
    if i > 0
      return i, j
    end
  end
  return i,-1
end
getIdx(pp::Symbol, sym::Symbol, i::Int=0)::Tuple{Int, Int} = pp==sym ? (abs(i)%100+1, div(abs(i)-100,100)) : (i-1, div(abs(i)-100,100))
function getIdx(pp::V, sym::Symbol, i::Int=0)::Tuple{Int, Int} where {V <: InferenceVariable}
  return getIdx(pp.dimtype, sym)
end

"""
    $SIGNATURES

Return whether `sym::Symbol` represents a variable vertex in the graph.
"""
isVariable(dfg::G, sym::Symbol) where G <: AbstractDFG = haskey(dfg.labelDict, sym)

"""
    $SIGNATURES

Return whether `sym::Symbol` represents a factor vertex in the graph.
"""
isFactor(dfg::G, sym::Symbol) where G <: AbstractDFG = hasFactor(dfg, sym)


# """
#     $SIGNATURES
#
# Return reference to a variable in `::FactorGraph` identified by `::Symbol`.
# """
# getVariable(fgl::FactorGraph, lbl::Symbol) = getVert(fgl, lbl, api=api)

# """
#     $SIGNATURES
#
# Return reference to the user factor in `::FactorGraph` identified by `::Symbol`.
# """
# getFactor(fvert::Graphs.ExVertex) = getData(fvert).fnc.usrfnc!
# getFactor(fgl::FactorGraph, lbl::Symbol, api::DataLayerAPI=dlapi) = getFactor(getVert(fgl, lbl, api=api, nt=:fct))

"""
    $SIGNATURES

Display and return to console the user factor identified by tag name.
"""
showFactor(fgl::G, fsym::Symbol) where G <: AbstractDFG = @show getFactor(fgl,fsym)

"""
   $SIGNATURES

Display the content of `VariableNodeData` to console for a given factor graph and variable tag`::Symbol`.
"""
function showVariable(fgl::G, vsym::Symbol) where G <: AbstractDFG
  vert = DFG.getVariable(fg, vsym)
  vnd = getData(vert)
  println("label: $(vert.label), exVertexId: $(vert.index)")
  println("tags: $( haskey(vert.attributes, string(:tags)) ? vert.attributes[string(:tags)] : string(:none))")
  println("size marginal samples $(size(getVal(vnd)))")
  println("kde bandwidths: $(getBW(vnd)[:,1])")
  println("kde mean: $(round.(getKDEMean(getKDE(vnd)),digits=4))")
  println("kde max: $(round.(getKDEMax(getKDE(vnd)),digits=4))")
  println()
  vnd
end

"""
    $SIGNATURES

Return `::Bool` on whether this variable has been marginalized.
"""
isMarginalized(vert::DFGVariable) = getData(vert).ismargin
isMarginalized(dfg::G, sym::Symbol; api::DataLayerAPI=localapi) where G <: AbstractDFG = isMarginalized(DFG.getVariable(fg, sym))



"""
    $SIGNATURES

Free all variables from marginalization.
"""
function unmarginalizeVariablesAll!(fgl::FactorGraph)
  fgl.solverParams.isfixedlag = false
  fgl.solverParams.qfl = 9999999999
  vsyms = union(ls(fgl)...)
  for sym in vsyms
    getData(fgl, sym).ismargin = false
  end
  nothing
end

"""
    $SIGNATURES

Free all variables from marginalization.

Related

unmarginalizeVariablesAll!
"""
unfreezeVariablesAll!(fgl::FactorGraph) = unmarginalizeVariablesAll!(fgl)

"""
    $SIGNATURES

Reset initialization flag on all variables in `::FactorGraphs`.

Notes
- Numerical values remain, but inference will overwrite since init flags are now `false`.
"""
function resetVariableAllInitializations!(fgl::FactorGraph)
  vsyms = union(ls(fgl)...)
  for sym in vsyms
    getData(fgl, sym).initialized = false
  end
  nothing
end



function convert(::Type{BallTreeDensity}, p::EasyMessage)
  AMP.manikde!(p.pts, p.bws, p.manifolds)
end

function convert(::Type{EasyMessage}, p::BallTreeDensity, manifolds::T) where {T <: Tuple}
  EasyMessage(getPoints(p), getBW(p)[:,1], manifolds)
end


#
