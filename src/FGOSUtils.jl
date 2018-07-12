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

function ls(fgl::FactorGraph)
  k = collect(keys(fgl.IDs))
  l = Int[]
  x = Int[]
  for kk in k
    kstr = string(kk)
    val = parse(Int,kstr[2:end]) # kk
    if kstr[1] == 'l'
      push!(l,val)
    elseif kstr[1] == 'x'
      push!(x,val)
    end
  end
  l = sort(l)
  x = sort(x)
  ll = Array{Symbol,1}(length(l))
  xx = Array{Symbol,1}(length(x))
  for i in 1:length(l)
    ll[i] = Symbol(string("l",l[i]))
  end
  for i in 1:length(x)
    xx[i] = Symbol(string("x",x[i]))
  end
  return xx,ll
end

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




#
