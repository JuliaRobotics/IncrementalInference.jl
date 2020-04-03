

function packmultihypo(fnc::CommonConvWrapper{T}) where {T<:FunctorInferenceType}
  @warn "packmultihypo is deprecated in favor of Vector only operations"
  fnc.hypotheses != nothing ? string(fnc.hypotheses) : ""
end
function parsemultihypostr(str::AS) where {AS <: AbstractString}
  @warn "parsemultihypostr is deprecated in favor of Vector only operations"
  mhcat=nothing
  if length(str) > 0
    mhcat = extractdistribution(str)
  end
  return mhcat
end


## packing converters-----------------------------------------------------------
# heavy use of multiple dispatch for converting between packed and original data types during DB usage

# function convert(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T}) where {P <: PackedInferenceType, T <: FunctorInferenceType}
#   # println("convert(::Type{PackedFunctionNodeData{$P}}, d::FunctionNodeData{$T})")
#   @error("convert GenericWrapParam is deprecated, use CommonConvWrapper instead.")
#   # mhstr = packmultihypo(d.fnc)
#   return PackedFunctionNodeData(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           string(d.frommodule), convert(P, d.fnc.usrfnc!), d.multihypo)
# end


function convert(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T}) where {P <: PackedInferenceType, T <: ConvolutionObject}
  # mhstr = packmultihypo(d.fnc)  # this is where certainhypo error occurs
  return PackedFunctionNodeData(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          string(d.frommodule), convert(P, d.fnc.usrfnc!),
          d.multihypo, d.fnc.certainhypo )  # extract two values from ccw for storage -- ccw thrown away
end



## unpack converters------------------------------------------------------------


function convert(
            ::Type{IncrementalInference.GenericFunctionNodeData{IncrementalInference.CommonConvWrapper{F},Symbol}},
            d::IncrementalInference.GenericFunctionNodeData{P,String} ) where {F <: FunctorInferenceType, P <: PackedInferenceType}
  #
  # TODO store threadmodel=MutliThreaded,SingleThreaded in persistence layer
  usrfnc = convert(F, d.fnc)
  # FIXME add proper nullhypo value
  mhcat, nh = parseusermultihypo(d.multihypo, 0.0)

  # TODO -- improve prepgenericconvolution for hypotheses and certainhypo field recovery when deserializing
  # reconstitute from stored data
  # FIXME, add threadmodel=threadmodel
  ccw = prepgenericconvolution(DFG.DFGVariable[], usrfnc, multihypo=mhcat, nullhypo=nh)
  ccw.certainhypo = d.certainhypo

  ret = FunctionNodeData{CommonConvWrapper{typeof(usrfnc)}}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs, Symbol(d.frommodule), ccw, d.multihypo, d.certainhypo)
  # error("what what $(ret.fnc.certainhypo)")
  return ret
end




function convert(::Type{PT}, ::T) where {PT <: PackedInferenceType, T <:FunctorInferenceType}
  getfield(T.name.module, Symbol("Packed$(T.name.name)"))
end
function convert(::Type{T}, ::PT) where {T <: FunctorInferenceType, PT <: PackedInferenceType}
  getfield(PT.name.module, Symbol(string(PT.name.name)[7:end]))
end


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
  cfnd = convert(PackedFunctionNodeData{usrtyp}, getSolverData(getFactor(fgl, fsym)) )
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
  neighborUserData = map(v->getSolverData(v).softtype, neighbors)

  # Rebuilding the CCW
  ccw_new = getDefaultFactorData(dfg, neighbors, getFactorType(factor), multihypo=getSolverData(factor).multihypo)
  setSolverData!(factor, ccw_new)

  #... Copying neighbor data into the factor?
  for i in 1:Threads.nthreads()
    ccw_new.fnc.cpt[i].factormetadata.variableuserdata = deepcopy(neighborUserData)
  end

  # Copying all other fields in the factor
  # TODO: Confirm whether we still need to do this?
  ## Rebuild getSolverData(fcnode).fncargvID, however, the list is order sensitive
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
  # fgs.g = Graphs.incdict(TreeClique,is_directed=false)

  # @showprogress 1 "Encoding variables..."
  for vsym in listVariables(fgl)
    # cpvert = deepcopy(  )
    var = getVariable(fgl, vsym)
    # addVariable!(fgs, cpvert)
  end

  # @showprogress 1 "Encoding factors..."
  for (fsym,fid) in fgs.fIDs
    data,ftyp = convert2packedfunctionnode(fgl, fsym)
    data = FunctionNodeData{ftyp}(Int[], false, false, Int[], m, gwpf)
    # newvert = TreeClique(fid,string(fsym))
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


# extend convenience function
function manikde!(pts::AbstractArray{Float64,2},
                  bws::Vector{Float64},
                  softtype::InferenceVariable  )::BallTreeDensity
  #
  manikde!(pts, bws, getManifolds(softtype))
end


function convert(::Type{Tuple{BallTreeDensity,Float64}},
                 p::TreeBelief )
  # @show size(p.val), size(p.bw), p.manifolds
  # (AMP.manikde!(p.val, p.bw[:,1], p.manifolds), p.inferdim)
  (convert(BallTreeDensity, p), p.inferdim)
end


function convert(::Type{TreeBelief},
                 bel::Tuple{BallTreeDensity,Float64},
                 manifolds::T) where {T <: Tuple}
  @error "Dont use this convert(::Type{TreeBelief}, bel::Tuple{BallTreeDensity,Float64}, manifolds) since it must assume ContinuousScalar softtype!!!"
  TreeBelief(getPoints(bel[1]), getBW(bel[1])[:,1:1], bel[2], ContinuousScalar(), manifolds)
end








#
