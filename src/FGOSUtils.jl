# Factor Graph OS type utilities
#  IIF methods should direclty detect extended types from user import
# of convert in their namespace


"""
    savefgjld(fgl::FactorGraph; file::AbstractString="tempfg.jld")

Save mostly complete Factor Graph type by converting complicated FunctionNodeData
types to 'Packed' types using user supplied converters. Ground truth can also be
saved and recovered by the associated loadjld(file="tempfg.jld") method.
"""
function savejld(fgl::FactorGraph;
      file::AbstractString="tempfg.jld",
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
    loadjld(file="tempfg.jld")

Opposite of savejld(fg, gt=gt, file="tempfg.jl") to load data from file. This function
uses the unpacking converters for converting all PackedInferenceType to FunctorInferenceType.
"""
function loadjld(;file::AbstractString="tempfg.jld")
  dd = JLD.load(file)
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
allnums(str::S) where {S <: AbstractString} = ismatch(Regex(string(["[0-9]" for j in 1:length(str)]...)), str)
# ismatch(r"_+|,+|-+", node_idx)

isnestednum(str::S; delim='_') where {S <: AbstractString} = ismatch(Regex("[0-9]+$(delim)[0-9]+"), str)

function sortnestedperm(strs::Vector{<:AbstractString}; delim='_')
  str12 = split.(strs, delim)
  sp1 = sortperm(parse.(Int,getindex.(str12,2)))
  sp2 = sortperm(parse.(Int,getindex.(str12,1)[sp1]))
  return sp1[sp2]
end


"""
    $(SIGNATURES)

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
lsf{T <: AbstractString}(fgl::FactorGraph, lbl::T) = lsf(fgl,Symbol(lbl))

function lsf{T <: FunctorInferenceType}(fgl::FactorGraph,
      mt::Type{T};
      api::DataLayerAPI=dlapi  )
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



#
