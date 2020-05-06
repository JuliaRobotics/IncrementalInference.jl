

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


function convert(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T}) where {P <: PackedInferenceType, T <: FactorOperationalMemory}
  # mhstr = packmultihypo(d.fnc)  # this is where certainhypo error occurs
  return PackedFunctionNodeData(d.eliminated, d.potentialused, d.edgeIDs,
          convert(P, d.fnc.usrfnc!),
          d.multihypo, d.fnc.certainhypo, d.solveInProgress)  # extract two values from ccw for storage -- ccw thrown away
end



## unpack converters------------------------------------------------------------


function convert(
            ::Type{IncrementalInference.GenericFunctionNodeData{IncrementalInference.CommonConvWrapper{F}}},
            d::IncrementalInference.GenericFunctionNodeData{P} ) where {F <: FunctorInferenceType, P <: PackedInferenceType}
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

  ret = FunctionNodeData{CommonConvWrapper{typeof(usrfnc)}}( d.eliminated, d.potentialused, d.edgeIDs, ccw, d.multihypo, d.certainhypo, d.solveInProgress)
  # error("what what $(ret.fnc.certainhypo)")
  return ret
end



# FIXME, move to DFG
function convert(::Type{PT}, ::T) where {PT <: PackedInferenceType, T <:FunctorInferenceType}
  getfield(T.name.module, Symbol("Packed$(T.name.name)"))
end
#FIXME This does not fit in with the convert funtion's expected behaviour
# it should rather be convert(::Type{T}, ::Type{PT}) or named something else
function convert(::Type{T}, ::PT) where {T <: FunctorInferenceType, PT <: PackedInferenceType}
  getfield(PT.name.module, Symbol(string(PT.name.name)[7:end]))
end


## Variables

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

  return factor
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





#
