
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


# typealias PackedFunctionNodeData{T <: PackedInferenceType} GenericFunctionNodeData{T, AbstractString}
PackedFunctionNodeData{T <: PackedInferenceType} = GenericFunctionNodeData{T, AbstractString}
PackedFunctionNodeData() = GenericFunctionNodeData{T, AbstractString}()
PackedFunctionNodeData(x1, x2, x3, x4, x5, x6) = GenericFunctionNodeData{T, AbstractString}(x1, x2, x3, x4, x5, x6)


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
#   return convert(P, d) #PackedVariableNodeData
# end
# function VNDdecoder(T::Type{VariableNodeData}, d::PackedVariableNodeData)
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


# heavy use of multiple dispatch for converting between packed and original data types during DB usage
function convert{T <: InferenceType, P <: PackedInferenceType}(::Type{FunctionNodeData{T}}, d::PackedFunctionNodeData{P})
  warn("P2P2BR should not be calling here")
  return FunctionNodeData{T}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          Symbol(d.frommodule), convert(T, d.fnc))
end
function convert{P <: PackedInferenceType, T <: InferenceType}(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T})
  return PackedFunctionNodeData{P}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          string(d.frommodule), convert(P, d.fnc))
end


# Functor version -- TODO, abstraction can be improved here
function convert(::Type{FunctionNodeData{GenericWrapParam{T}}},
            d::PackedFunctionNodeData{P} ) where {T <: FunctorInferenceType, P <: PackedInferenceType}
  #
  # @show "convert", T, P
  # @show typeof(d.fnc)
  # info("calling convert($(T), $(d.fnc))")
  usrfnc = convert(T, d.fnc)
  # @show typeof(usrfnc)
  gwpf = prepgenericwrapper(Graphs.ExVertex[], usrfnc, getSample)
  return FunctionNodeData{GenericWrapParam{typeof(usrfnc)}}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          Symbol(d.frommodule), gwpf) #{T}
end
function convert{P <: PackedInferenceType, T <: FunctorInferenceType}(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T})
  return PackedFunctionNodeData{P}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          string(d.frommodule), convert(P, d.fnc.usrfnc!))
end

# function FNDencode{T <: FunctorInferenceType, P <: PackedInferenceType}(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T})
#   return convert(PackedFunctionNodeData{P}, d) #PackedFunctionNodeData{P}
# end
# function FNDdecode{T <: FunctorInferenceType, P <: PackedInferenceType}(::Type{FunctionNodeData{T}}, d::PackedFunctionNodeData{P})
#   return convert(FunctionNodeData{T}, d) #FunctionNodeData{T}
# end
#
# function FNDencode{T <: InferenceType, P <: PackedInferenceType}(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T})
#   return convert(PackedFunctionNodeData{P}, d) #PackedFunctionNodeData{P}
# end
# function FNDdecode{T <: InferenceType, P <: PackedInferenceType}(::Type{FunctionNodeData{T}}, d::PackedFunctionNodeData{P})
#   return convert(FunctionNodeData{T}, d) #FunctionNodeData{T}
# end


# Compare FunctionNodeData
function compare{T,S}(a::GenericFunctionNodeData{T,S},b::GenericFunctionNodeData{T,S})
  # TODO -- beef up this comparison to include the gwp
  TP = true
  TP = TP && a.fncargvID == b.fncargvID
  TP = TP && a.eliminated == b.eliminated
  TP = TP && a.potentialused == b.potentialused
  TP = TP && a.edgeIDs == b.edgeIDs
  TP = TP && a.frommodule == b.frommodule
  TP = TP && typeof(a.fnc) == typeof(b.fnc)
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
function encodePackedType(topackdata::FunctionNodeData{T}) where {T}
  println("IncrementalInference.encodePackedType(::PackedFunctionNodeData{T}): Unknown packed type encoding T=$(T) of $(topackdata)")
  # @show T, typeof(topackdata)
  fnctype = getfnctype(topackdata)
  # @show sfnctype = split(string(fnctype),'.')
  fnc = getfield(getmodule(fnctype), Symbol("Packed$(getname(fnctype))"))
  convert(PackedFunctionNodeData{fnc}, topackdata)
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
function decodePackedType(packeddata::GenericFunctionNodeData{PT,S}, typestring::String) where {PT, S <: AbstractString}
  # @show typestring
  # @show typeof(packeddata), PT
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
