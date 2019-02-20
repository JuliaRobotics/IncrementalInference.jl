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
    compareAll(Al[i], Bl[i], show=show, skip=skip)
  end
  return true
end

function compareAll(Al::T,
                    Bl::T;
                    show::Bool=true,
                    skip::Vector{Symbol}=Symbol[]  )::Bool where {T <: DataType}
  #
  return true
end

function compareAll(Al::T, Bl::T; show::Bool=true, skip::Vector{Symbol}=Symbol[])::Bool where T
  TP = compareFields(Al, Bl, show=show, skip=skip)
  for field in fieldnames(T)
    if field in skip
      continue
    end
    Ad = eval(:($Al.$field))
    Bd = eval(:($Bl.$field))
    TP = TP && compareAll(Ad, Bd, show=show, skip=skip)
  end
  return TP
end


function compare(p1::BallTreeDensity, p2::BallTreeDensity)::Bool
  return compareAll(p1.bt,p2.bt, skip=[:calcStatsHandle; :data]) &&
         compareAll(p1,p2, skip=[:calcStatsHandle; :bt])
end

"""
    $(SIGNATURES)

Save mostly complete Factor Graph type by converting complicated FunctionNodeData
types to 'Packed' types using user supplied converters. Ground truth can also be
saved and recovered by the associated loadjld(file="tempfg.jld2") method.
"""
function savejld(fgl::FactorGraph;
      file::AbstractString="tempfg.jld2",
      groundtruth=nothing)
  fgs = encodefg(fgl)
  if groundtruth == nothing
    @save file fgs
  else
    @save file fgs groundtruth
  end
  return file
end



"""
    $(SIGNATURES)

Opposite of savejld(fg, gt=gt, file="tempfg.jl") to load data from file. This function
uses the unpacking converters for converting all PackedInferenceType to FunctorInferenceType.
"""
function loadjld(;file::AbstractString="tempfg.jld2")
  dd = FileIO.load(file)
  # fgs = jldopen(file,"r") do file
  #   read(file, "fgs")
  # end
  fgs, gt = nothing, nothing
  if haskey(dd, "fgs")
    fgs = dd["fgs"]
  else
    error("No factor graph (fgs) data found in this file, only found $(keys(dd))")
  end
  if haskey(dd, "groundtruth")
    gt = dd["groundtruth"]
    println("Also found ground truth data")
  end
  fgd = decodefg(fgs)
  return fgd, gt
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

Return all elements `ls(fg)` as tuples, or nodes connected to the a specific element, eg. `ls(fg, :x1)
"""
function ls(fgl::FactorGraph, lbl::Symbol; api::DataLayerAPI=dlapi, ring::Int=1)
  # TODO ring functionality must still be implemented
  lsa = Symbol[]
  # v = nothing
  if haskey(fgl.IDs, lbl)
    id = fgl.IDs[lbl]
  else
    return lsa
  end
  # this is unnecessary
  v = getVert(fgl,id, api=api)
  for outn in api.outneighbors(fgl, v)
    # if outn.attributes["ready"] = 1 && outn.attributes["backendset"]=1
      push!(lsa, Symbol(outn.label))
    # end
  end
  return lsa
end
ls(fgl::FactorGraph, lbl::T) where {T <: AbstractString} = ls(fgl, Symbol(lbl))

"""
    $(SIGNATURES)

Experimental union of elements version of ls(::FactorGraph, ::Symbol).  Not mean't to replace broadcasting `ls.(fg, [:x1;:x2])`
"""
function ls(fgl::FactorGraph,
            lbls::Vector{Symbol};
            api::DataLayerAPI=dlapi,
            ring::Int=1)
  union(ls.(fgl, lbls, ring=ring, api=api)[:]...)
end

"""
    $(SIGNATURES)

List the nodes in a factor graph.

# Examples
```julia-repl
ls(fg)
```
"""
function ls(fgl::FactorGraph; key1='x', key2='l')
  k = collect(keys(fgl.IDs))
  x, l = String[], String[]
  xval, lval = Int[], Int[]
  xstr, lstr = String[], String[]
  xvalnested, lvalnested = String[], String[]
  xstrnested, lstrnested = String[], String[]
  canparse1, canparse2 = true,true
  nestedparse1, nestedparse2 = true, true
  idx = 0
  for id in k
    idx += 1
    idstr = string(id)
    # val = parse(Int,kstr[2:end]) # TODO: handle non-int labels
    node_idx = idstr[2:end]
    canparse = allnums(node_idx)
    nested = isnestednum(node_idx)
    if idstr[1] == key1
      keystr = string(key1,node_idx)
      if canparse
        push!(xstr, keystr)
        push!(xval, parse(Int, node_idx))
      elseif nested
        push!(xvalnested, node_idx)
        push!(xstrnested, string(node_idx))
      else
        push!(x,keystr)
      end
    elseif idstr[1] == key2
      keystr = string(key2,node_idx)
      if canparse
        push!(lstr, keystr)
        push!(lval, parse(Int, node_idx))
      elseif nested
        push!(lstrnested, keystr)
        push!(lvalnested, string(node_idx))
      else
        push!(l,string(key2,node_idx))
      end
    end
  end
  x1 = xstr[sortperm(xval)]
  x2 = xstrnested[sortnestedperm(xvalnested)]
  x = [x1; x2; sort(x)]

  l1 = lstr[sortperm(lval)]
  l2 = lstrnested[sortnestedperm(lvalnested)]
  l = [l1; l2; sort(l)]

  xx = Symbol.(x)
  ll = Symbol.(l)
  return xx, ll #return poses, landmarks
end

"""
    $(SIGNATURES)

List factors in a factor graph.

# Examples
```julia-repl
lsf(fg)
```
"""
function lsf(fgl::FactorGraph, lbl::Symbol; api::DataLayerAPI=dlapi)
  lsa = Symbol[]
  # v = Union{}
  if haskey(fgl.fIDs, lbl)
    id = fgl.fIDs[lbl]
  else
    return lsa
  end
  v = getVert(fgl, id, api=api) # fgl.g.vertices[id] #fgl.f[id]
  for outn in api.outneighbors(fgl, v) # out_neighbors(v, fgl.g)
    push!(lsa, Symbol(outn.label))
  end
  return lsa
end

"""
    $(SIGNATURES)

List factors in a factor graph.

# Examples
```julia-repl
lsf(fg)
```
"""
lsf(fgl::FactorGraph, lbl::T) where {T <: AbstractString} = lsf(fgl,Symbol(lbl))

"""
    $(SIGNATURES)

List factors in a factor graph.

# Examples
```julia-repl
lsf(fg)
```
"""
function lsf(fgl::FactorGraph,
      mt::Type{T};
      api::DataLayerAPI=dlapi  ) where {T <: FunctorInferenceType}
  #
  syms = Symbol[]
  for (fsym,fid) in fgl.fIDs
    if typeof(getfnctype(fgl, fid, api=api))==T
      push!(syms, fsym)
    end
  end
  return syms
end

function ls2(fgl::FactorGraph, vsym::Symbol)
  xxf = ls(fgl, vsym)
  xlxl = Symbol[]
  for xf in xxf
    xx = lsf(fgl,xf)
    xlxl = union(xlxl, xx)
  end
  xlxl = setdiff(xlxl, [vsym])
  return xlxl
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


hasOrphans(fg) = sum(length.(ls.(fg, [ls(fg)[1];ls(fg)[2]])) .== 0) > 0

"""
    landmarks(fgl::FactorGraph, vsym::Symbol)

Return Vector{Symbol} of landmarks attached to vertex vsym in fgl.
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

Return the user factor identified by tag name.
"""
getFactor(fgl::FactorGraph, fsym::Symbol) = getData(fgl, fsym, nt=:fct).fnc.usrfnc!

"""
    $SIGNATURES

Display and return to console the user factor identified by tag name.
"""
showFactor(fgl::FactorGraph, fsym::Symbol; api::DataLayerAPI=dlapi) = @show getFactor(fgl,fsym)

"""
   $SIGNATURES

Display the content of `VariableNodeData` to console for a given factor graph and variable tag`::Symbol`. 
"""
function showVariable(fgl::FactorGraph, vsym::Symbol; api::DataLayerAPI=dlapi)
  vert = getVert(fg, vsym, api=api)
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

#
