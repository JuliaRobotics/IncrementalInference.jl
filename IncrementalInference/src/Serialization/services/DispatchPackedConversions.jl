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
  dfg::AbstractDFG,
  factor::FactorCompute,
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
  #   mergeState!(factor, new_solverData)
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
