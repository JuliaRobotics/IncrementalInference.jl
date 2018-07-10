
# Creation of packed types should be automated with macros

mutable struct PackedVariableNodeData
  vecinitval::Array{Float64,1}
  diminitval::Int
  vecinitstdev::Array{Float64,1}
  diminitdev::Int
  vecval::Array{Float64,1}
  dimval::Int
  vecbw::Array{Float64,1}
  dimbw::Int
  BayesNetOutVertIDs::Array{Int,1}
  dimIDs::Array{Int,1}
  dims::Int
  eliminated::Bool
  BayesNetVertID::Int
  separator::Array{Int,1}
  # groundtruth::VoidUnion{ Dict{ Tuple{Symbol, Vector{Float64}} } }
  softtype::String
  initialized::Bool
  PackedVariableNodeData() = new()
  PackedVariableNodeData(x1::Vector{Float64},
                         x2::Int,
                         x3::Vector{Float64},
                         x4::Int,
                         x5::Vector{Float64},
                         x6::Int,
                         x7::Vector{Float64},
                         x8::Int,
                         x9::Vector{Int},
                         x10::Vector{Int},
                         x11::Int,
                         x12::Bool,
                         x13::Int,
                         x14::Vector{Int},
                         x15::String,
                         x16::Bool) = new(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16)
end



const FunctionNodeData{T <: Union{InferenceType, FunctorInferenceType}} = GenericFunctionNodeData{T, Symbol}
# FunctionNodeData() = GenericFunctionNodeData{T, Symbol}()
FunctionNodeData(x1, x2, x3, x4, x5::Symbol, x6::T, x7::String="") where {T <: FunctorInferenceType}= GenericFunctionNodeData{T, Symbol}(x1, x2, x3, x4, x5, x6, x7)

# typealias PackedFunctionNodeData{T <: PackedInferenceType} GenericFunctionNodeData{T, AbstractString}
const PackedFunctionNodeData{T <: PackedInferenceType} = GenericFunctionNodeData{T, <: AbstractString}
# PackedFunctionNodeData{T}() where T = GenericFunctionNodeData{T, AbstractString}()
# PackedFunctionNodeData(x1, x2, x3, x4, x5::S, x6::T) where {T <: PackedInferenceType, S <: AbstractString} = GenericFunctionNodeData{T, AbstractString}(x1, x2, x3, x4, x5, x6)
PackedFunctionNodeData(x1, x2, x3, x4, x5::S, x6::T, x7::String="") where {T <: PackedInferenceType, S <: AbstractString} = GenericFunctionNodeData(x1, x2, x3, x4, x5, x6, x7)



function normalfromstring(str::AS) where {AS <: AbstractString}
  meanstr = match(r"μ=[+-]?([0-9]*[.])?[0-9]+", str).match
  mean = split(meanstr, '=')[2]
  sigmastr = match(r"σ=[+-]?([0-9]*[.])?[0-9]+", str).match
  sigma = split(sigmastr, '=')[2]
  Normal{Float64}(parse(Float64,mean), parse(Float64,sigma))
end

function mvnormalfromstring(str::AS) where {AS <: AbstractString}
  means = split(split(split(str, 'μ')[2],']')[1],'[')[end]
  mean = Float64[]
  for ms in split(means, ',')
    push!(mean, parse(Float64, ms))
  end
  sigs = split(split(split(str, 'Σ')[2],']')[1],'[')[end]
  sig = Float64[]
  for ms in split(sigs, ';')
    for m in split(ms, ' ')
      length(m) > 0 ? push!(sig, parse(Float64, m)) : nothing
    end
  end
  len = length(mean)
  sigm = reshape(sig, len,len)
  MvNormal(mean, sigm)
end

function categoricalfromstring(str::AS)::Distributions.Categorical where {AS <: AbstractString}
  # pstr = match(r"p=\[", str).match
  psubs = split(str, '=')[end]
  psubs = split(psubs, '[')[end]
  psubsub = split(psubs, ']')[1]
  pw = split(psubsub, ',')
  return Categorical(parse.(Float64, pw))
end

function extractdistribution(str::AS)::Union{Void, Distributions.Distribution} where {AS <: AbstractString}
  if str == ""
    return nothing
  elseif (ismatch(r"Normal", str) && !ismatch(r"FullNormal", str))
    return normalfromstring(str)
  elseif ismatch(r"FullNormal", str)
    return mvnormalfromstring(str)
  elseif ismatch(r"Categorical", str)
    return categoricalfromstring(str)
  else
    error("Don't know how to extract distrubtion from str=$(str)")
  end
end


function convert(::Type{PackedVariableNodeData}, d::VariableNodeData)
  return PackedVariableNodeData(d.initval[:],size(d.initval,1),
                              d.initstdev[:],size(d.initstdev,1),
                              d.val[:],size(d.val,1),
                              d.bw[:], size(d.bw,1),
                              d.BayesNetOutVertIDs,
                              d.dimIDs, d.dims, d.eliminated,
                              d.BayesNetVertID, d.separator,
                              string(d.softtype), d.initialized)
end
function convert(::Type{VariableNodeData}, d::PackedVariableNodeData)

  r1 = d.diminitval
  c1 = r1 > 0 ? floor(Int,length(d.vecinitval)/r1) : 0
  M1 = reshape(d.vecinitval,r1,c1)

  r2 = d.diminitdev
  c2 = r2 > 0 ? floor(Int,length(d.vecinitstdev)/r2) : 0
  M2 = reshape(d.vecinitstdev,r2,c2)

  r3 = d.dimval
  c3 = r3 > 0 ? floor(Int,length(d.vecval)/r3) : 0
  M3 = reshape(d.vecval,r3,c3)

  r4 = d.dimbw
  c4 = r4 > 0 ? floor(Int,length(d.vecbw)/r4) : 0
  M4 = reshape(d.vecbw,r4,c4)

  # TODO -- allow out of module type allocation (future feature, not currently in use)
  st = IncrementalInference.ContinuousMultivariate # eval(parse(d.softtype))

  return VariableNodeData(M1,M2,M3,M4, d.BayesNetOutVertIDs,
    d.dimIDs, d.dims, d.eliminated, d.BayesNetVertID, d.separator,
    nothing, st, d.initialized )
end
# function VNDencoder(P::Type{PackedVariableNodeData}, d::VariableNodeData)
#   warn("VNDencoder deprecated, use the convert functions through dispatch instead, P=$(P).")
#   return convert(P, d) #PackedVariableNodeData
# end
# function VNDdecoder(T::Type{VariableNodeData}, d::PackedVariableNodeData)
#   warn("VNDdecoder deprecated, use the convert functions through dispatch instead, T=$(T).")
#   return convert(T, d) #VariableNodeData
# end


function compare(a::VariableNodeData,b::VariableNodeData)
    TP = true
    TP = TP && a.initval == b.initval
    TP = TP && a.initstdev == b.initstdev
    TP = TP && a.val == b.val
    TP = TP && a.bw == b.bw
    TP = TP && a.BayesNetOutVertIDs == b.BayesNetOutVertIDs
    TP = TP && a.dimIDs == b.dimIDs
    TP = TP && a.dims == b.dims
    TP = TP && a.eliminated == b.eliminated
    TP = TP && a.BayesNetVertID == b.BayesNetVertID
    TP = TP && a.separator == b.separator
    return TP
end

function ==(a::VariableNodeData,b::VariableNodeData, nt::Symbol=:var)
  return IncrementalInference.compare(a,b)
end

# TODO change to rather use the extractdistribution functions that already exists
function packmultihypo(fnc::GenericWrapParam{T}) where {T<:FunctorInferenceType}
  # fnc.hypotheses != nothing ? "$(fnc.hypotheses.p)" : ""
  fnc.hypotheses != nothing ? string(fnc.hypotheses) : ""
end
function parsemultihypostr(str::AS) where {AS <: AbstractString}
  # mhverts = Symbol[]
  mhcat=nothing
  if length(str) > 0
    mhcat = extractdistribution(str)
    # # ss = split(str, ';')
    # # s1 = strip.(split(split(split(ss[1],']')[1],'[')[end], ','))
    # s2 = strip.(split(split(split(str,']')[1],'[')[end], ','))
    # # mhverts = Symbol.(String[s1[i][1]==':'? s1[i][2:end] : s1[i] for i in 1:length(s1)])
    # mhcat = Distributions.Categorical(parse.(Float64, s2))
  end
  return mhcat
end

# heavy use of multiple dispatch for converting between packed and original data types during DB usage


# TODO simplify and reduce to single pack and unpack converter
## packing converters-----------------------------------------------------------

function convert(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T}) where {P <: PackedInferenceType, T <: InferenceType}
    # println("convert(::Type{PackedFunctionNodeData{$P}}, d::FunctionNodeData{$T})")
  error("convert(::Type{PackedFunctionNodeData{P}} : this convert function should probably not be used. Use FunctorInferenceType instead of InferenceType.")
  mhstr = packmultihypo(d.fnc)
  return PackedFunctionNodeData{P}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          string(d.frommodule), convert(P, d.fnc), mhstr) # TODO -- should use convert(P, d.fnc.usrfnc!) instead
end
function convert(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T}) where {P <: PackedInferenceType, T <: FunctorInferenceType}
  # println("convert(::Type{PackedFunctionNodeData{$P}}, d::FunctionNodeData{$T})")
  mhstr = packmultihypo(d.fnc)
  return PackedFunctionNodeData(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          string(d.frommodule), convert(P, d.fnc.usrfnc!), mhstr)
end
# function convert(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T}) where {P, T <: PackedInferenceType}
#   # TODO weird to have functions call this converter, but okay for now...
#   println("convert(::Type{PackedFunctionNodeData{$P}}, d::FunctionNodeData{$T})")
#   return PackedFunctionNodeData(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           string(d.frommodule), d.fnc.usrfnc!)
# end




## unpack converters------------------------------------------------------------

# Functor version -- TODO, abstraction can be improved here
function convert(::Type{FunctionNodeData{T}}, d::PackedFunctionNodeData{P}) where {T <: InferenceType, P <: PackedInferenceType}
  error("Old unpacking converter from DB to Graphs.jl")
  return FunctionNodeData{T}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          Symbol(d.frommodule), convert(T, d.fnc))
end
# function convert(::Type{FunctionNodeData{GenericWrapParam{T}}},
#             d::PackedFunctionNodeData{P} ) where {T <: FunctorInferenceType, P <: PackedInferenceType}
#   #
#   warn("Unpacking Option 1, F=$(F), P=$(P)")
#   usrfnc = convert(T, d.fnc)
#   gwpf = prepgenericwrapper(Graphs.ExVertex[], usrfnc, getSample)
#   return FunctionNodeData{GenericWrapParam{typeof(usrfnc)}}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           Symbol(d.frommodule), gwpf) #{T}
# end
function convert(
            ::Type{IncrementalInference.GenericFunctionNodeData{IncrementalInference.GenericWrapParam{F},Symbol}},
            d::IncrementalInference.GenericFunctionNodeData{P,String} ) where {F <: FunctorInferenceType, P <: PackedInferenceType}
  #
  # warn("Unpacking Option 2, F=$(F), P=$(P)")
  usrfnc = convert(F, d.fnc)
  # @show d.multihypo
  mhcat = parsemultihypostr(d.multihypo)
  # @show typeof(mhcat)
  gwpf = prepgenericwrapper(Graphs.ExVertex[], usrfnc, getSample, multihypo=mhcat)
  return FunctionNodeData{GenericWrapParam{typeof(usrfnc)}}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          Symbol(d.frommodule), gwpf)
end






function FNDencode{T <: FunctorInferenceType, P <: PackedInferenceType}(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T})
  warn("FNDencode deprecated, use the convert functions through dispatch instead, PackedFunctionNodeData{P=$(P)}.")
  return convert(PackedFunctionNodeData{P}, d) #PackedFunctionNodeData{P}
end
function FNDdecode{T <: FunctorInferenceType, P <: PackedInferenceType}(::Type{FunctionNodeData{T}}, d::PackedFunctionNodeData{P})
  warn("FNDdecode deprecated, use the convert functions through dispatch instead, FunctionNodeData{T=$(T)}.")
  return convert(FunctionNodeData{T}, d) #FunctionNodeData{T}
end

function FNDencode{T <: InferenceType, P <: PackedInferenceType}(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T})
  warn("FNDencode deprecated, use the convert functions through dispatch instead, PackedFunctionNodeData{P=$(P)}.")
  return convert(PackedFunctionNodeData{P}, d) #PackedFunctionNodeData{P}
end
function FNDdecode{T <: InferenceType, P <: PackedInferenceType}(::Type{FunctionNodeData{T}}, d::PackedFunctionNodeData{P})
  warn("FNDdecode deprecated, use the convert functions through dispatch instead, FunctionNodeData{T=$(T)}.")
  return convert(FunctionNodeData{T}, d) #FunctionNodeData{T}
end


# Compare FunctionNodeData
function compare(a::GenericFunctionNodeData{T1,S},b::GenericFunctionNodeData{T2,S}) where {T1, T2, S}
  # TODO -- beef up this comparison to include the gwp
  TP = true
  TP = TP && a.fncargvID == b.fncargvID
  TP = TP && a.eliminated == b.eliminated
  TP = TP && a.potentialused == b.potentialused
  TP = TP && a.edgeIDs == b.edgeIDs
  TP = TP && a.frommodule == b.frommodule
  # TP = TP && typeof(a.fnc) == typeof(b.fnc)
  return TP
end



function convert{PT <: PackedInferenceType, T <:FunctorInferenceType}(::Type{PT}, ::T)
  getfield(T.name.module, Symbol("Packed$(T.name.name)"))
end
function convert{T <: FunctorInferenceType, PT <: PackedInferenceType}(::Type{T}, ::PT)
  getfield(PT.name.module, Symbol(string(PT.name.name)[7:end]))
end



# function encodePackedType(topackdata::T) where T
#   println("IncrementalInference.encodePackedType(data) new dispatch conversion to packed type development:")
#   @show fnc = getfield(T.name.module, Symbol("Packed$(T.name.name)"))
#   return convert(fnc, topackdata)
#   # error("IncrementalInference.encodePackedType(::$(fnc)) unknown format")
#   # data = nothing
#   # if fnc == PackedFunctionNodeData
#   #   @show usrtyp = convert(PackedInferenceType, fnc)
#   #   @show data = convert(IncrementalInference.PackedFunctionNodeData{usrtyp}, topackdata )
#   # else
#   #
#   # end
#   # return data
# end

# function encodePackedType(topackdata::GenericWrapParam{T}) where {T}
#   @show T
#   error("encodePackedType")
# end
function getmodule(t::T) where T
  T.name.module
end
function getname(t::T) where T
  T.name.name
end
function encodePackedType(topackdata::VariableNodeData)
  # error("IncrementalInference.encodePackedType(::VariableNodeData): Unknown packed type encoding of $(topackdata)")
  convert(IncrementalInference.PackedVariableNodeData, topackdata)
end
function encodePackedType(topackdata::GenericFunctionNodeData{T, Symbol}) where {T <: FunctorInferenceType}
  # println("IncrementalInference.encodePackedType(::GenericFunctionNodeData{T,Symbol}): Unknown packed type encoding T=$(T) of $(topackdata)")
  fnctype = getfnctype(topackdata)
  fnc = getfield(getmodule(fnctype), Symbol("Packed$(getname(fnctype))"))
  convert(PackedFunctionNodeData{fnc}, topackdata)
end
function encodePackedType(topackdata::GenericFunctionNodeData{T, <:AbstractString}) where {T <: PackedInferenceType}
  error("IncrementalInference.encodePackedType(::FunctionNodeData{T, <:AbstractString}): Unknown packed type encoding T=$(T) of $(topackdata)")
  # @show T, typeof(topackdata)
  # warn("Yes, its packed!")
  # fnctype = getfnctype(topackdata)
  # @show fnc = getfield(getmodule(fnctype), Symbol("Packed$(getname(fnctype))"))
  # convert(PackedFunctionNodeData{T}, topackdata)
  topackdata
end

"""
    $(SIGNATURES)

Encode complicated function node type to related 'Packed<type>' format assuming a user supplied convert function .
"""
function convert2packedfunctionnode(fgl::FactorGraph,
      fsym::Symbol,
      api::DataLayerAPI=localapi  )
  #
  fid = fgl.fIDs[fsym]
  fnc = getfnctype(fgl, fid)
  usrtyp = convert(PackedInferenceType, fnc)
  cfnd = convert(PackedFunctionNodeData{usrtyp},getData(fgl, fid, api=api) )
  return cfnd, usrtyp
end

# function encodePackedType(topackdata::T) where {T <: FunctorInferenceType}
#   error("IncrementalInference.encodePackedType(<: FunctorInferenceType): Unknown packed type encoding of $(fnc)")
# end


#


# function decodePackedType(packeddata::PT, typestring::String) where {PT}
#   @show typeof(packeddata), PT
#   @show tse = split(typestring,'.')[end]  # TODO -- make more robust by checking the modules as well
#   @show fulltype = getfield(PT.name.module, Symbol(tse[7:end]))
#   # @show packedtype = eval(parse(typestring))
#   return convert(fulltype, packeddata)
#
#   # if packedtype == PackedVariableNodeData
#   #   return convert(VariableNodeData, packeddata)
#   # else
#   #   error("IncrementalInference.decodePackedType doesnt know how to handle $(packedtype) yet")
#   # end
# end
function decodePackedType(packeddata::PackedVariableNodeData, typestring::String)
  # error("IncrementalInference.encodePackedType(::VariableNodeData): Unknown packed type encoding of $(topackdata)")
  convert(IncrementalInference.VariableNodeData, packeddata)
end
function decodePackedType(packeddata::GenericFunctionNodeData{PT,<:AbstractString}, typestring::String) where {PT}
  # warn("decodePackedType($(typeof(packeddata)),$(typestring)) is happening with PT=$(PT) and ")
  functype = getfield(PT.name.module, Symbol(string(PT.name.name)[7:end]))
  fulltype = FunctionNodeData{GenericWrapParam{functype}}
  convert(fulltype, packeddata)
end

"""
    $(SIGNATURES)

Make a full memory copy of the graph and encode all complicated function node
types with assumed to exist convert to 'Packed<type>' formats. Same converters
as used for database persistence storage with CloudGraphs.jl.
"""
function encodefg(fgl::FactorGraph;
      api::DataLayerAPI=localapi  )
  #
  fgs = deepcopy(fgl)
  fgs.cg = nothing
  fgs.registeredModuleFunctions = nothing
  fgs.g = Graphs.incdict(Graphs.ExVertex,is_directed=false)
  @showprogress 1 "Encoding variables..." for (vsym,vid) in fgl.IDs
    cpvert = deepcopy(getVert(fgl, vid, api=api))
    api.addvertex!(fgs, cpvert) #, labels=vnlbls)  # currently losing labels
  end

  @showprogress 1 "Encoding factors..." for (fsym,fid) in fgs.fIDs
    data,ftyp = convert2packedfunctionnode(fgl, fsym)
    # data = FunctionNodeData{ftyp}(Int[], false, false, Int[], m, gwpf)
    newvert = ExVertex(fid,string(fsym))
    for (key,val) in getVert(fgl,fid,api=api).attributes
      newvert.attributes[key] = val
    end
    setData!(newvert, data)
    api.addvertex!(fgs, newvert)
  end
  fgs.g.inclist = typeof(fgl.g.inclist)()

  # iterated over all edges
  @showprogress 1 "Encoding edges..." for (eid, edges) in fgl.g.inclist
    fgs.g.inclist[eid] = Vector{typeof(edges[1])}()
    for ed in edges
      newed = Graphs.Edge(ed.index,
          fgs.g.vertices[ed.source.index],
          fgs.g.vertices[ed.target.index]  )
      push!(fgs.g.inclist[eid], newed)
    end
  end

  return fgs
end


"""
    convertfrompackedfunctionnode(fgl, fsym)

If you get unknown type conversion error when loading a .jld, while using your own
FunctorInferenceTypes, you should:
Copy these functions below, and overload in your package with local extented
FunctorInferenceType definitions.
See RoME/src/fgos.jl for example.
"""
function convertfrompackedfunctionnode(fgl::FactorGraph,
      fsym::Symbol,
      api::DataLayerAPI=localapi  )
  #
  fid = fgl.fIDs[fsym]
  fnc = getData(fgl, fid).fnc #getfnctype(fgl, fid)
  usrtyp = convert(FunctorInferenceType, fnc)
  data = getData(fgl, fid, api=api)
  newtype = FunctionNodeData{GenericWrapParam{usrtyp}}
  cfnd = convert(newtype, data)
  return cfnd, usrtyp
end

"""
    decodefg(fgs::FactorGraph)

Unpack PackedFunctionNodeData formats back to regular FunctonNodeData.
"""
function decodefg(fgs::FactorGraph; api::DataLayerAPI=localapi)
  fgu = deepcopy(fgs)
  fgu.cg = nothing
  fgu.registeredModuleFunctions = nothing
  fgu.g = Graphs.incdict(Graphs.ExVertex,is_directed=false)
  @showprogress 1 "Decoding variables..." for (vsym,vid) in fgs.IDs
    cpvert = deepcopy(getVert(fgs, vid, api=api))
    api.addvertex!(fgu, cpvert) #, labels=vnlbls)  # currently losing labels
  end

  @showprogress 1 "Decoding factors..." for (fsym,fid) in fgu.fIDs
    data,ftyp = convertfrompackedfunctionnode(fgs, fsym)
    # data = FunctionNodeData{ftyp}(Int[], false, false, Int[], m, gwpf)
    newvert = ExVertex(fid,string(fsym))
    for (key,val) in getVert(fgs,fid,api=api).attributes
      newvert.attributes[key] = val
    end
    setData!(newvert, data)
    api.addvertex!(fgu, newvert)
  end
  fgu.g.inclist = typeof(fgs.g.inclist)()

  # iterated over all edges
  @showprogress 1 "Decoding edges..." for (eid, edges) in fgs.g.inclist
    fgu.g.inclist[eid] = Vector{typeof(edges[1])}()
    for ed in edges
      newed = Graphs.Edge(ed.index,
          fgu.g.vertices[ed.source.index],
          fgu.g.vertices[ed.target.index]  )
      push!(fgu.g.inclist[eid], newed)
    end
  end

  return fgu
end





# MethodError: Cannot `convert` an object of type
#
# IncrementalInference.GenericFunctionNodeData{RoME.PackedPriorPose2,String} to an object of type
# IncrementalInference.GenericFunctionNodeData{IncrementalInference.GenericWrapParam{RoME.PriorPose2},Symbol}
#
# This may have arisen from a call to the constructor
# IncrementalInference.GenericFunctionNodeData{IncrementalInference.GenericWrapParam{RoME.PriorPose2},Symbol}(...),
# since type constructors fall back to convert methods.
# decodePackedType(::IncrementalInference.GenericFunctionNodeData{RoME.PackedPriorPose2,String}, ::String) at DispatchPackedConversions.jl:269
# unpackNeoNodeData2UsrType(::CloudGraphs.CloudGraph, ::Neo4j.Node) at CloudGraphs.jl:160
# neoNode2CloudVertex(::CloudGraphs.CloudGraph, ::Neo4j.Node) at CloudGraphs.jl:171
# get_vertex(::CloudGraphs.CloudGraph, ::Int64, ::Bool) at CloudGraphs.jl:242
# makeAddCloudEdge!(::IncrementalInference.FactorGraph, ::Graphs.ExVertex, ::Graphs.ExVertex) at CloudGraphIntegration.jl:231
# #addFactor!#27(::Int64, ::IncrementalInference.DataLayerAPI, ::Array{String,1}, ::Int64, ::Bool, ::Function, ::IncrementalInference.FactorGraph, ::Array{Graphs.ExVertex,1}, ::RoME.PriorPose2) at FactorGraph01.jl:519
# (::IncrementalInference.#kw##addFactor!)(::Array{Any,1}, ::IncrementalInference.#addFactor!, ::IncrementalInference.FactorGraph, ::Array{Graphs.ExVertex,1}, ::RoME.PriorPose2) at <missing>:0
# #addFactor!#28(::Int64, ::IncrementalInference.DataLayerAPI, ::Array{String,1}, ::Int64, ::Bool, ::Function, ::IncrementalInference.FactorGraph, ::Array{Symbol,1}, ::RoME.PriorPose2) at FactorGraph01.jl:546
# addFactor!(::IncrementalInference.FactorGraph, ::Array{Symbol,1}, ::RoME.PriorPose2) at FactorGraph01.jl:542
# include_string(::String, ::String) at loading.jl:522
# include_string(::String, ::String, ::Int64) at eval.jl:30
# include_string(::Module, ::String, ::String, ::Int64, ::Vararg{Int64,N} where N) at eval.jl:34
# (::Atom.##102#107{String,Int64,String})() at eval.jl:82
# withpath(::Atom.##102#107{String,Int64,String}, ::String) at utils.jl:30
# withpath(::Function, ::String) at eval.jl:38
# hideprompt(::Atom.##101#106{String,Int64,String}) at repl.jl:67
# macro expansion at eval.jl:80 [inlined]
# (::Atom.##100#105{Dict{String,Any}})() at task.jl:80














#
