# functions relating to parametric solutions of a single factor that is likely in need of consolidation

"""
    $SIGNATURES

Helper function to propagate a parametric estimate along a factor chain.  
This function takes and returns variable values as coordinates.

Notes
- Not used during MM-iSAM inference.
- Expected uses are for user analysis of factors and estimates.
- real-time dead reckoning chain prediction.
- Parametric binary factor utility function, used by DRT

DevNotes
- TODO ensure type stability, likely returning types `Any` at this time.
- TODO MeanMaxPPE currently stored as coordinates, complicating fast calculation.

Related: [`getMeasurementParametric`](@ref), [`approxConvBelief`](@ref), [`MutablePose2Pose2Gaussian`](@ref)
"""
function solveFactorParametric(
  dfg::AbstractDFG,
  fct::DFGFactor,
  # currval::P1,
  srcsym_vals::AbstractVector{Pair{Symbol, P}},
  trgsym::Symbol;
  solveKey::Symbol = :default,
  evaltmpkw...,
) where {P}
  #

  varLbls = getVariableOrder(fct)
  varTypes = tuple((getVariableType.(dfg, varLbls))...)
  sfidx = findfirst(varLbls .== trgsym)

  # get the measurement point
  fctTyp = getFactorType(fct)
  # this is definitely in coordinates, see JuliaRobotics/RoME.jl#465
  mea, _ = getMeasurementParametric(fctTyp)
  # must change measT to be a tangent vector
  M = getManifold(fctTyp)
  e0 = getPointIdentity(M)
  mea_ = hat(M, e0, mea)
  measT = [mea_]

  # get variable points
  function _getParametric(vari::DFGVariable, key = :default)
    # hasp = haskey(getPPEDict(vari), key)
    # FIXME use PPE via Manifold points currently in coordinates
    # hasp ? getPPE(vari, key).suggested : calcMean(getBelief(vari, key))
    pt = calcMean(getBelief(vari, key))

    return collect(getCoordinates(getVariableType(vari), pt))
  end

  # overwrite specific src values from user
  coordVals = _getParametric.(getVariable.(dfg, varLbls), solveKey)

  for (srcsym, currval) in srcsym_vals
    coordVals[findfirst(varLbls .== srcsym)] = currval
  end
  crds = tuple(coordVals...)

  pts = tuple(map(t -> getPoint(t...), zip(varTypes, crds))...)

  # do the calculation to find solvefor index using the factor, as manifold point
  pt = _evalFactorTemporary!(
    fctTyp,
    varTypes,
    sfidx,
    measT,
    pts;
    solveKey,
    evaltmpkw...,
  )[1]

  return getCoordinates(getVariableType(dfg, trgsym), pt)
end
