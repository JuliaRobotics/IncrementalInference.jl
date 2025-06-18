
## packing converters-----------------------------------------------------------
# heavy use of multiple dispatch for converting between packed and original data types during DB usage

function convert(
  ::Type{PackedFunctionNodeData{P}},
  d::FunctionNodeData{T},
) where {P <: AbstractPackedFactor, T <: FactorSolverCache}
  error("TODO remove. PackedFunctionNodeData is obsolete")
  return PackedFunctionNodeData(
    d.eliminated,
    d.potentialused,
    d.edgeIDs,
    convert(P, _getCCW(d).usrfnc!),
    d.multihypo,
    _getCCW(d).hyporecipe.certainhypo,
    d.nullhypo,
    d.solveInProgress,
    d.inflation,
  )  # extract two values from ccw for storage -- ccw thrown away
end

## unpack converters------------------------------------------------------------

# see #1424
#TODO Consolidate: this looks alot like `getDefaultFactorData`
function DFG.reconstFactorData(
  dfg::AbstractDFG,
  varOrder::AbstractVector{Symbol},
  ::Type{<:GenericFunctionNodeData{<:CommonConvWrapper{F}}},
  packed::GenericFunctionNodeData{<:AbstractPackedFactor},
) where {F <: AbstractFactor}

  error("TODO remove. Obsolete: use `DFG.rebuildFactorCache!` and getDefaultFactorData instead.")
  #
  # TODO store threadmodel=MutliThreaded,SingleThreaded in persistence layer
  usrfnc = convert(F, packed.fnc)
  multihypo, nullhypo = parseusermultihypo(packed.multihypo, packed.nullhypo)

  # IIF #1424
  vars = map(f -> getVariable(dfg, f), varOrder)
  userCache = preambleCache(dfg, vars, usrfnc)

  # TODO -- improve _createCCW for hypotheses and certainhypo field recovery when deserializing
  # reconstitute from stored data
  # FIXME, add threadmodel=threadmodel
  # FIXME https://github.com/JuliaRobotics/DistributedFactorGraphs.jl/issues/590#issuecomment-776838053
  # FIXME dont know what manifolds to use in ccw
  ccw = _createCCW(
    vars,
    usrfnc;
    multihypo,
    nullhypo,
    certainhypo = packed.certainhypo,
    inflation = packed.inflation,
    userCache,
    attemptGradients = getSolverParams(dfg).attemptGradients,
    # Block recursion if NoSolverParams or if set to not attempt gradients.
    _blockRecursion=
      getSolverParams(dfg) isa NoSolverParams || 
      !getSolverParams(dfg).attemptGradients,
  )
  #

  # CommonConvWrapper{typeof(usrfnc)}
  ret = FunctionNodeData{typeof(ccw)}(
    packed.eliminated,
    packed.potentialused,
    packed.edgeIDs,
    ccw,
    packed.multihypo,
    packed.certainhypo,
    packed.nullhypo,
    packed.solveInProgress,
    packed.inflation,
  )
  #
  return ret
end

##

"""
    $(SIGNATURES)

After deserializing a factor using decodePackedType, use this to
completely rebuild the factor's CCW and user data.

Notes:
- This function is likely to be used for cache heavy factors, e.g. `ObjectAffordanceSubcloud`.

Dev Notes:
- TODO: We should only really do this in-memory if we can by without it (review this).
- TODO: needs testing
"""
function DFG.rebuildFactorCache!(
  dfg::AbstractDFG{SolverParams},
  factor::DFGFactor,
  neighbors = map(vId -> getVariable(dfg, vId), listNeighbors(dfg, factor));
  _blockRecursionGradients::Bool=false
)
  #
  # Set up the neighbor data

  # Rebuilding the CCW
  state = DFG.getFactorState(factor)
  state, solvercache = getDefaultFactorData(
    dfg,
    neighbors,
    DFG.getObservation(factor);
    multihypo = state.multihypo,
    nullhypo = state.nullhypo,
    # special inflation override
    inflation = state.inflation,
    eliminated = state.eliminated,
    potentialused = state.potentialused,
    solveInProgress = state.solveInProgress,
    _blockRecursion=_blockRecursionGradients
  )
  #
  DFG.setCache!(factor, solvercache)
  return factor

  # factor_ = if typeof(solvercache) != typeof(DFG.getCache(factor)) 
  #   # must change the type of factor solver data FND{CCW{...}}
  #   # create a new factor
  #   factor__ = FactorCompute(
  #     getLabel(factor),
  #     Tuple(getVariableOrder(factor)),
  #     DFG.getObservation(factor),
  #     state,
  #     solvercache;
  #     timestamp = getTimestamp(factor),
  #     nstime = factor.nstime,
  #     tags = getTags(factor),
  #     solvable = getSolvable(factor),
  #   )
  #   #

  #   # replace old factor in dfg with a new one
  #   deleteFactor!(dfg, factor; suppressGetFactor = true)
  #   addFactor!(dfg, factor__)

  #   factor__
  # else
  #   setSolverData!(factor, new_solverData)
  #   DFG.setCache!(factor, solvercache)
  #   # We're not updating here because we don't want
  #   # to solve cloud in loop, we want to make sure this flow works:
  #   # Pull big cloud graph into local -> solve local -> push back into cloud.
  #   # updateFactor!(dfg, factor)
  #   factor
  # end

  #... Copying neighbor data into the factor?
  # JT TODO it looks like this is already updated in getDefaultFactorData -> _createCCW
  # factormetadata.variableuserdata is deprecated, remove when removing deprecation
  # for i in 1:Threads.nthreads()
  #   ccw_new.fnc.cpt[i].factormetadata.variableuserdata = deepcopy(neighborUserData)
  # end

  # return factor_
end

## =================================================================
## TODO Can code below be deprecated?
## =================================================================

function convert(::Type{Tuple{ManifoldKernelDensity, Float64}}, p::TreeBelief)
  # 
  return (convert(ManifoldKernelDensity, p), p.infoPerCoord)
end

#
