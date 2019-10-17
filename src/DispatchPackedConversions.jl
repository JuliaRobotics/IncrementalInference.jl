

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
  @error("convert GenericWrapParam is deprecated, use CommonConvWrapper instead.")
  mhstr = packmultihypo(d.fnc)
  return PackedFunctionNodeData(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          string(d.frommodule), convert(P, d.fnc.usrfnc!), mhstr)
end


function convert(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T}) where {P <: PackedInferenceType, T <: ConvolutionObject}
  mhstr = packmultihypo(d.fnc)  # this is where certainhypo error occurs
  return PackedFunctionNodeData(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          string(d.frommodule), convert(P, d.fnc.usrfnc!),
          mhstr, d.fnc.certainhypo )  # extract two values from ccw for storage -- ccw thrown away
end



## unpack converters------------------------------------------------------------


function convert(
            ::Type{IncrementalInference.GenericFunctionNodeData{IncrementalInference.CommonConvWrapper{F},Symbol}},
            d::IncrementalInference.GenericFunctionNodeData{P,String} ) where {F <: FunctorInferenceType, P <: PackedInferenceType}
  #
  # TODO store threadmodel=MutliThreaded,SingleThreaded in persistence layer
  usrfnc = convert(F, d.fnc)
  mhcat = parsemultihypostr(d.multihypo)

  # TODO -- improve prepgenericconvolution for hypotheses and certainhypo field recovery when deserializing
  # reconstitute from stored data
  ccw = prepgenericconvolution(DFG.DFGVariable[], usrfnc, multihypo=mhcat)
  ccw.certainhypo = d.certainhypo

  ret = FunctionNodeData{CommonConvWrapper{typeof(usrfnc)}}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          Symbol(d.frommodule), ccw)
  # error("what what $(ret.fnc.certainhypo)")
  return ret
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

# function getmodule(t::T) where T
#   T.name.module
# end
# function getname(t::T) where T
#   T.name.name
# end
# function getpackedtype(typestring::AS) where {AS <: AbstractString}
#   # println("Caesar.getpackedtype($(typestring))")
#   eval(Meta.parse(typestring))() # TODO consider caching or better
# end
# function encodePackedType(topackdata::VariableNodeData)
#   @warn "'encodePackedType' Deprecated..."
#   # error("IncrementalInference.encodePackedType(::VariableNodeData): Unknown packed type encoding of $(topackdata)")
#   pack(nothing, topackdata)
# end
# function encodePackedType(topackdata::GenericFunctionNodeData{CommonConvWrapper{T}, Symbol}) where {T <: FunctorInferenceType}
#   # println("IncrementalInference.encodePackedType(::GenericFunctionNodeData{T,Symbol}): Unknown packed type encoding T=$(T) of $(topackdata)")
#   fnctype = getfnctype(topackdata)
#   fnc = getfield(getmodule(fnctype), Symbol("Packed$(getname(fnctype))"))
#   convert(PackedFunctionNodeData{fnc}, topackdata)
# end
# function encodePackedType(topackdata::GenericFunctionNodeData{T, <:AbstractString}) where {T <: PackedInferenceType}
#   error("IncrementalInference.encodePackedType(::FunctionNodeData{T, <:AbstractString}): Unknown packed type encoding T=$(T) of $(topackdata)")
#   # @show T, typeof(topackdata)
#   # @warn "Yes, its packed!"
#   # fnctype = getfnctype(topackdata)
#   # @show fnc = getfield(getmodule(fnctype), Symbol("Packed$(getname(fnctype))"))
#   # convert(PackedFunctionNodeData{T}, topackdata)
#   topackdata
# end

"""
    $(SIGNATURES)

Encode complicated function node type to related 'Packed<type>' format assuming a user supplied convert function .
"""
function convert2packedfunctionnode(fgl::G,
                                    fsym::Symbol ) where G <: AbstractDFG
  #
  # fid = fgl.fIDs[fsym]
  fnc = getfnctype(fgl, fsym)
  usrtyp = convert(PackedInferenceType, fnc)
  cfnd = convert(PackedFunctionNodeData{usrtyp}, solverData(getFactor(fgl, fsym)) )
  return cfnd, usrtyp
end



# Variables

function decodePackedType(dfg::G, packeddata::PackedVariableNodeData) where G <: AbstractDFG
  @warn "Deprecated?"
  convert(IncrementalInference.VariableNodeData, packeddata)
end
# Factors
function decodePackedType(dfg::G, packeddata::GenericFunctionNodeData{PT,<:AbstractString}) where {PT, G <: AbstractDFG}
  usrtyp = convert(FunctorInferenceType, packeddata.fnc)
  fulltype = FunctionNodeData{CommonConvWrapper{usrtyp}}
  factor = convert(fulltype, packeddata)
  return factor
end

"""
  $(SIGNATURES)
After deserializing a factor using decodePackedType, use this to
completely rebuild the factor's CCW and user data.
"""
function rebuildFactorMetadata!(dfg::G, factor::DFGFactor)::DFGFactor where G <: AbstractDFG
  # Set up the neighbor data
  neighbors = map(vId->getVariable(dfg, vId), getNeighbors(dfg, factor))
  neighborUserData = map(v->solverData(v).softtype, neighbors)

  # Rebuilding the CCW
  setDefaultFactorNode!(dfg, factor, neighbors, factor.data.fnc.usrfnc!)

  #... Copying neighbor data into the factor?
  ccw_new = solverData(factor)
  for i in 1:Threads.nthreads()
    ccw_new.fnc.cpt[i].factormetadata.variableuserdata = deepcopy(neighborUserData)
  end

  # Copying all other fields in the factor
  # TODO: Confirm whether we still need to do this?
  ## Rebuild solverData(fcnode).fncargvID, however, the list is order sensitive
  # out_neighbors does not gaurantee ordering -- i.e. why is it not being saved
  # for field in fieldnames(typeof(ccw_jld))
  #   if field != :fnc
  #     setfield!(ccw_new, field, getfield(ccw_jld, field))
  #   end
  # end

  return factor
end



"""
    $(SIGNATURES)

Make a full memory copy of the graph and encode all composite function node
types -- assuming that convert methods for 'Packed<type>' formats exist.  The same converters
are used for database persistence with CloudGraphs.jl.
"""
function encodefg(fgl::G ) where G <: AbstractDFG
  #
  fgs = deepcopy(fgl)
  # fgs.g = Graphs.incdict(Graphs.ExVertex,is_directed=false)

  # @showprogress 1 "Encoding variables..."
  for vsym in getVariableIds(fgl)
    # cpvert = deepcopy(  )
    var = getVariable(fgl, vsym)
    # addVariable!(fgs, cpvert)
  end

  # @showprogress 1 "Encoding factors..."
  for (fsym,fid) in fgs.fIDs
    data,ftyp = convert2packedfunctionnode(fgl, fsym)
    data = FunctionNodeData{ftyp}(Int[], false, false, Int[], m, gwpf)
    # newvert = ExVertex(fid,string(fsym))
    # for (key,val) in getVert(fgl,fid,api=api).attributes
    #   newvert.attributes[key] = val
    # end
    ## losing fgl.fncargvID before setdata
    # setData!(newvert, data)
    # api.addvertex!(fgs, newvert)
  end
  fgs.g.inclist = typeof(fgl.g.inclist)()

  # iterated over all edges
  # @showprogress 1 "Encoding edges..."
  for (eid, edges) in fgl.g.inclist
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

# import IncrementalInference: decodefg, loadjld

veeCategorical(val::Categorical) = val.p
veeCategorical(val::Union{Nothing, Vector{Float64}}) = val


















#
