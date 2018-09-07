
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
  ismargin::Bool
  dontmargin::Bool
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
                         x16::Bool,
                         x17::Bool,
                         x18::Bool ) = new(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18)
end



# # TODO stop-gap string storage of Distrubtion types, should be upgraded to more efficient storage
# function normalfromstring(str::AS) where {AS <: AbstractString}
#   meanstr = match(r"μ=[+-]?([0-9]*[.])?[0-9]+", str).match
#   mean = split(meanstr, '=')[2]
#   sigmastr = match(r"σ=[+-]?([0-9]*[.])?[0-9]+", str).match
#   sigma = split(sigmastr, '=')[2]
#   Normal{Float64}(parse(Float64,mean), parse(Float64,sigma))
# end
#
# function mvnormalfromstring(str::AS) where {AS <: AbstractString}
#   means = split(split(split(str, 'μ')[2],']')[1],'[')[end]
#   mean = Float64[]
#   for ms in split(means, ',')
#     push!(mean, parse(Float64, ms))
#   end
#   sigs = split(split(split(str, 'Σ')[2],']')[1],'[')[end]
#   sig = Float64[]
#   for ms in split(sigs, ';')
#     for m in split(ms, ' ')
#       length(m) > 0 ? push!(sig, parse(Float64, m)) : nothing
#     end
#   end
#   len = length(mean)
#   sigm = reshape(sig, len,len)
#   MvNormal(mean, sigm)
# end
#
# function categoricalfromstring(str::AS)::Distributions.Categorical where {AS <: AbstractString}
#   # pstr = match(r"p=\[", str).match
#   psubs = split(str, '=')[end]
#   psubs = split(psubs, '[')[end]
#   psubsub = split(psubs, ']')[1]
#   pw = split(psubsub, ',')
#   p = parse.(Float64, pw)
#   return Categorical(p ./ sum(p))
# end
#
# function extractdistribution(str::AS)::Union{Void, Distributions.Distribution} where {AS <: AbstractString}
#   # TODO improve use of multidispatch and packing of Distribution types
#   if str == ""
#     return nothing
#   elseif (ismatch(r"Normal", str) && !ismatch(r"FullNormal", str))
#     return normalfromstring(str)
#   elseif ismatch(r"FullNormal", str)
#     return mvnormalfromstring(str)
#   elseif ismatch(r"Categorical", str)
#     return categoricalfromstring(str)
#   elseif ismatch(r"KDE:", str)
#     return convert(KDE.BallTreeDensity, str)
#   elseif ismatch(r"AliasingScalarSampler", str)
#     return convert(AliasingScalarSampler, str)
#   else
#     error("Don't know how to extract distribution from str=$(str)")
#   end
# end


function convert(::Type{PackedVariableNodeData}, d::VariableNodeData)
  return PackedVariableNodeData(d.initval[:],size(d.initval,1),
                                d.initstdev[:],size(d.initstdev,1),
                                d.val[:],size(d.val,1),
                                d.bw[:], size(d.bw,1),
                                d.BayesNetOutVertIDs,
                                d.dimIDs, d.dims, d.eliminated,
                                d.BayesNetVertID, d.separator,
                                string(d.softtype), d.initialized, d.ismargin, d.dontmargin)
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
    nothing, st, d.initialized, d.ismargin, d.dontmargin )
end



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
    TP = TP && a.ismargin == b.ismargin
    return TP
end

function ==(a::VariableNodeData,b::VariableNodeData, nt::Symbol=:var)
  return IncrementalInference.compare(a,b)
end


function packmultihypo(fnc::CommonConvWrapper{T}) where {T<:FunctorInferenceType}
  fnc.hypotheses != nothing ? string(fnc.hypotheses) : ""
end
function parsemultihypostr(str::AS) where {AS <: AbstractString}
  mhcat=nothing
  if length(str) > 0
    mhcat = extractdistribution(str)
  end
  return mhcat
end


## packing converters-----------------------------------------------------------
# heavy use of multiple dispatch for converting between packed and original data types during DB usage

function convert(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T}) where {P <: PackedInferenceType, T <: FunctorInferenceType}
  # println("convert(::Type{PackedFunctionNodeData{$P}}, d::FunctionNodeData{$T})")
  warn("convert GenericWrapParam is deprecated, use CommonConvWrapper instead.")
  mhstr = packmultihypo(d.fnc)
  return PackedFunctionNodeData(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          string(d.frommodule), convert(P, d.fnc.usrfnc!), mhstr)
end
function convert(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T}) where {P <: PackedInferenceType, T <: ConvolutionObject}
  # println("convert(::Type{PackedFunctionNodeData{$P}}, d::FunctionNodeData{$T})")
  mhstr = packmultihypo(d.fnc)
  return PackedFunctionNodeData(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          string(d.frommodule), convert(P, d.fnc.usrfnc!), mhstr)
end



## unpack converters------------------------------------------------------------


function convert(
            ::Type{IncrementalInference.GenericFunctionNodeData{IncrementalInference.CommonConvWrapper{F},Symbol}},
            d::IncrementalInference.GenericFunctionNodeData{P,String} ) where {F <: FunctorInferenceType, P <: PackedInferenceType}
  #
  # warn("Unpacking Option 2, F=$(F), P=$(P)")
  usrfnc = convert(F, d.fnc)
  # @show d.multihypo
  mhcat = parsemultihypostr(d.multihypo)
  # TODO store threadmodel=MutliThreaded,SingleThreaded in persistence layer
  ccw = prepgenericconvolution(Graphs.ExVertex[], usrfnc, multihypo=mhcat)
  return FunctionNodeData{CommonConvWrapper{typeof(usrfnc)}}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          Symbol(d.frommodule), ccw)
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



function convert(::Type{PT}, ::T) where {PT <: PackedInferenceType, T <:FunctorInferenceType}
  getfield(T.name.module, Symbol("Packed$(T.name.name)"))
end
function convert(::Type{T}, ::PT) where {T <: FunctorInferenceType, PT <: PackedInferenceType}
  getfield(PT.name.module, Symbol(string(PT.name.name)[7:end]))
end



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
  cfnd = convert(PackedFunctionNodeData{usrtyp}, getData(fgl, fid, api=api) )
  return cfnd, usrtyp
end



function decodePackedType(packeddata::PackedVariableNodeData,
                          typestring::String )
  #
  # error("IncrementalInference.encodePackedType(::VariableNodeData): Unknown packed type encoding of $(topackdata)")
  convert(IncrementalInference.VariableNodeData, packeddata)
end
function decodePackedType(packeddata::GenericFunctionNodeData{PT,<:AbstractString}, typestring::String) where {PT}
  # warn("decodePackedType($(typeof(packeddata)),$(typestring)) is happening with PT=$(PT) and ")
  functype = getfield(PT.name.module, Symbol(string(PT.name.name)[7:end]))
  @assert convert(FunctorInferenceType, PT) == functype # NEW METHOD: TODO -- use as primary once known to work consistently
  fulltype = FunctionNodeData{CommonConvWrapper{functype}}
  convert(fulltype, packeddata)
end

"""
    $(SIGNATURES)

Make a full memory copy of the graph and encode all composite function node
types -- assuming that convert methods for 'Packed<type>' formats exist.  The same converters
are used for database persistence with CloudGraphs.jl.
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
    ## losing fgl.fncargvID before setdata
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
    $(SIGNATURES)

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
  newtype = FunctionNodeData{CommonConvWrapper{usrtyp}}
  cfnd = convert(newtype, data)
  return cfnd, usrtyp
end

# import IncrementalInference: decodefg, loadjld

veeCategorical(val::Categorical) = val.p
veeCategorical(val::Union{Void, Vector{Float64}}) = val


"""
    $(SIGNATURES)

Unpack PackedFunctionNodeData formats back to regular FunctonNodeData.
"""
function decodefg(fgs::FactorGraph; api::DataLayerAPI=localapi)
  fgu = deepcopy(fgs)
  fgu.cg = nothing # will be deprecated or replaced
  fgu.registeredModuleFunctions = nothing # obsolete
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

  # rebuild factormetadata
  @showprogress 1 "Rebuilding factor metadata..." for (fsym,fid) in fgu.fIDs
    varuserdata = []
    fcnode = getVert(fgu, fsym, nt=:fnc)
    # ccw = getData(fcnode)
    ccw_jld = deepcopy(getData(fcnode))
    allnei = Graphs.ExVertex[]
    for nei in out_neighbors(fcnode, fgu.g)
        push!(allnei, nei)
        data = IncrementalInference.getData(nei)
        push!(varuserdata, data.softtype)
    end
    setDefaultFactorNode!(fgu, fcnode, allnei, ccw_jld.fnc.usrfnc!, threadmodel=ccw_jld.fnc.threadmodel, multihypo=veeCategorical(ccw_jld.fnc.hypotheses))
    ccw_new = IncrementalInference.getData(fcnode)
    for i in 1:Threads.nthreads()
      ccw_new.fnc.cpt[i].factormetadata.variableuserdata = deepcopy(varuserdata)
    end
    ## Rebuild getData(fcnode).fncargvID, however, the list is order sensitive
    # out_neighbors does not gaurantee ordering -- i.e. why is it not being saved
    for field in fieldnames(ccw_jld)
      if field != :fnc
        setfield!(ccw_new, field, getfield(ccw_jld, field))
      end
    end
  end
  return fgu
end


















#
