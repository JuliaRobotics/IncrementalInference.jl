

function packmultihypo(fnc::CommonConvWrapper{T}) where {T<:AbstractFactor}
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

# FIXME see #1424
function convert(
            ::Type{GenericFunctionNodeData{CommonConvWrapper{F}}},
            packed::GenericFunctionNodeData{P} ) where {F <: AbstractFactor, P <: PackedInferenceType}
  #
  # TODO store threadmodel=MutliThreaded,SingleThreaded in persistence layer
  usrfnc = convert(F, packed.fnc)
  mhcat, nh = parseusermultihypo(packed.multihypo, packed.nullhypo)

  # TODO -- improve _prepCCW for hypotheses and certainhypo field recovery when deserializing
  # reconstitute from stored data
  # FIXME, add threadmodel=threadmodel
  # FIXME https://github.com/JuliaRobotics/DistributedFactorGraphs.jl/issues/590#issuecomment-776838053
  # FIXME dont know what manifolds to use in ccw
  ccw = _prepCCW(DFG.DFGVariable[], usrfnc, multihypo=mhcat, nullhypo=nh, inflation=packed.inflation)
  ccw.certainhypo = packed.certainhypo

  # CommonConvWrapper{typeof(usrfnc)}
  ret = FunctionNodeData{typeof(ccw)}(packed.eliminated, packed.potentialused, packed.edgeIDs, ccw,
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
Dev Notes:
- TODO: We should only really do this in-memory if we can by without it (review this).
"""
function rebuildFactorMetadata!(dfg::AbstractDFG{SolverParams}, 
                                factor::DFGFactor,
                                neighbors = map(vId->getVariable(dfg, vId), getNeighbors(dfg, factor)) )
  #
  # Set up the neighbor data

  # Rebuilding the CCW
  fsd = getSolverData(factor)
  fnd_new = getDefaultFactorData( dfg, 
                                  neighbors, 
                                  getFactorType(factor), 
                                  multihypo=fsd.multihypo,
                                  nullhypo=fsd.nullhypo,
                                  # special inflation override 
                                  inflation=fsd.inflation,
                                  eliminated=fsd.eliminated,
                                  potentialused=fsd.potentialused,
                                  edgeIDs=fsd.edgeIDs,
                                  solveInProgress=fsd.solveInProgress)
  #
  
  factor_ = if typeof(fnd_new) != typeof(getSolverData(factor))
    # must change the type of factor solver data FND{CCW{...}}
    # create a new factor
    factor__ = DFGFactor(getLabel(factor),
                        getTimestamp(factor),
                        factor.nstime,
                        getTags(factor),
                        fnd_new,
                        getSolvable(factor),
                        Tuple(getVariableOrder(factor)))
    #

    # replace old factor in dfg with a new one
    deleteFactor!(dfg, factor, suppressGetFactor=true)
    addFactor!(dfg, factor__)

    factor__
  else
    setSolverData!(factor, fnd_new)
    # We're not updating here because we don't want
    # to solve cloud in loop, we want to make sure this flow works:
    # Pull big cloud graph into local -> solve local -> push back into cloud.
    # updateFactor!(dfg, factor)
    factor
  end

  #... Copying neighbor data into the factor?
  # JT TODO it looks like this is already updated in getDefaultFactorData -> _prepCCW
  # factormetadata.variableuserdata is deprecated, remove when removing deprecation
  # for i in 1:Threads.nthreads()
  #   ccw_new.fnc.cpt[i].factormetadata.variableuserdata = deepcopy(neighborUserData)
  # end

  return factor_
end


# import IncrementalInference: decodefg, loadjld

veeCategorical(val::Categorical) = val.p
veeCategorical(val::Union{Nothing, Vector{Float64}}) = val


function convert( ::Type{Tuple{ManifoldKernelDensity,Float64}},
                  p::TreeBelief )
  # 
  (convert(ManifoldKernelDensity, p), p.infoPerCoord)
end





#
