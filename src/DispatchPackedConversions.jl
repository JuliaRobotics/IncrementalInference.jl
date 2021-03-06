

function packmultihypo(fnc::CommonConvWrapper{T}) where {T<:FunctorInferenceType}
  @warn "packmultihypo is deprecated in favor of Vector only operations"
  fnc.hypotheses !== nothing ? string(fnc.hypotheses) : ""
end
function parsemultihypostr(str::AS) where {AS <: AbstractString}
  @warn "parsemultihypostr is deprecated in favor of Vector only operations"
  mhcat=nothing
  if length(str) > 0
    mhcat = convert(SamplableBelief, str)
  end
  return mhcat
end


## packing converters-----------------------------------------------------------
# heavy use of multiple dispatch for converting between packed and original data types during DB usage


function convert(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T}) where {P <: PackedInferenceType, T <: FactorOperationalMemory}
  return PackedFunctionNodeData(d.eliminated, d.potentialused, d.edgeIDs,
                                convert(P, _getCCW(d).usrfnc!),
                                d.multihypo, _getCCW(d).certainhypo, d.nullhypo, 
                                d.solveInProgress, d.inflation)  # extract two values from ccw for storage -- ccw thrown away
end



## unpack converters------------------------------------------------------------


function convert(
            ::Type{GenericFunctionNodeData{CommonConvWrapper{F}}},
            packed::GenericFunctionNodeData{P} ) where {F <: FunctorInferenceType, P <: PackedInferenceType}
  #
  # TODO store threadmodel=MutliThreaded,SingleThreaded in persistence layer
  usrfnc = convert(F, packed.fnc)
  mhcat, nh = parseusermultihypo(packed.multihypo, packed.nullhypo)

  # TODO -- improve prepgenericconvolution for hypotheses and certainhypo field recovery when deserializing
  # reconstitute from stored data
  # FIXME, add threadmodel=threadmodel
  # FIXME https://github.com/JuliaRobotics/DistributedFactorGraphs.jl/issues/590#issuecomment-776838053
  ccw = prepgenericconvolution(DFG.DFGVariable[], usrfnc, multihypo=mhcat, nullhypo=nh, inflation=packed.inflation)
  ccw.certainhypo = packed.certainhypo

  ret = FunctionNodeData{CommonConvWrapper{typeof(usrfnc)}}(packed.eliminated, packed.potentialused, packed.edgeIDs, ccw,
                                                            packed.multihypo, packed.certainhypo, packed.nullhypo, 
                                                            packed.solveInProgress, packed.inflation )
  #
  return ret
end



## Variables

"""
  $(SIGNATURES)
After deserializing a factor using decodePackedType, use this to
completely rebuild the factor's CCW and user data.
"""
function rebuildFactorMetadata!(dfg::AbstractDFG{SolverParams}, factor::DFGFactor)
  # Set up the neighbor data
  neighbors = map(vId->getVariable(dfg, vId), getNeighbors(dfg, factor))
  neighborUserData = map(v->getVariableType(v), neighbors)

  # Rebuilding the CCW
  fsd = getSolverData(factor)
  ccw_new = getDefaultFactorData( dfg, neighbors, getFactorType(factor), 
                                  multihypo=fsd.multihypo,
                                  nullhypo=fsd.nullhypo,
                                  # special inflation override 
                                  inflation=fsd.inflation,
                                  eliminated=fsd.eliminated,
                                  potentialused=fsd.potentialused,
                                  edgeIDs=fsd.edgeIDs,
                                  solveInProgress=fsd.solveInProgress)
  setSolverData!(factor, ccw_new)

  #... Copying neighbor data into the factor?
  # JT TODO it looks like this is already updated in getDefaultFactorData -> prepgenericconvolution
  # factormetadata.variableuserdata is deprecated, remove when removing deprecation
  # for i in 1:Threads.nthreads()
  #   ccw_new.fnc.cpt[i].factormetadata.variableuserdata = deepcopy(neighborUserData)
  # end

  return factor
end


# import IncrementalInference: decodefg, loadjld

veeCategorical(val::Categorical) = val.p
veeCategorical(val::Union{Nothing, Vector{Float64}}) = val


function convert( ::Type{Tuple{BallTreeDensity,Float64}},
                  p::TreeBelief )
  # @show size(p.val), size(p.bw), p.manifolds
  # (AMP.manikde!(p.val, p.bw[:,1], p.manifolds), p.inferdim)
  (convert(BallTreeDensity, p), p.inferdim)
end





#
