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
  neighbors = collect(getVariable.(dfg, getVariableOrder(factor)));
  _blockRecursionGradients::Bool=false
)
  prepareFactorCache!(dfg, factor, neighbors; _blockRecursion = _blockRecursionGradients)
  return nothing
end

## =================================================================
## TODO Can code below be deprecated?
## =================================================================

function convert(::Type{Tuple{ManifoldKernelDensity, Float64}}, p::TreeBelief)
  # 
  return (convert(ManifoldKernelDensity, p), p.infoPerCoord)
end

#
