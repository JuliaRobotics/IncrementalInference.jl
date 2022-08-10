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
- TODO consolidate with [`solveConditionalsParametric`](@ref)

See also: [`getMeasurementParametric`](@ref), [`solveConditionalsParametric`](@ref), [`approxConvBelief`](@ref), [`MutablePose2Pose2Gaussian`](@ref)
"""
function solveFactorParameteric(dfg::AbstractDFG,
                                fct::DFGFactor,
                                # currval::P1,
                                srcsym_vals::AbstractVector{Pair{Symbol, P}},
                                trgsym::Symbol,
                                solveKey::Symbol=:default;
                                evaltmpkw...  ) where P
  #

  varLbls = getVariableOrder(fct)
  varTypes = tuple((getVariableType.(dfg, varLbls))...)
  sfidx = findfirst( varLbls .== trgsym )

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
  function _getParametric(vari::DFGVariable, key=:default)
    # hasp = haskey(getPPEDict(vari), key)
    # FIXME use PPE via Manifold points currently in coordinates
    # hasp ? getPPE(vari, key).suggested : calcMean(getBelief(vari, key))
    pt = calcMean(getBelief(vari, key))

    getCoordinates(getVariableType(vari),pt)
  end

  # overwrite specific src values from user
  coordVals = _getParametric.(getVariable.(dfg, varLbls), solveKey)
  for (srcsym, currval) in srcsym_vals
    coordVals[findfirst(varLbls .== srcsym)] = currval
  end
  crds = tuple(coordVals...)
  
  pts = tuple(map(t->getPoint(t...), zip(varTypes,crds))...)

  # do the calculation to find solvefor index using the factor, as manifold point
  pt = _evalFactorTemporary!( fctTyp, varTypes, sfidx, measT, pts; solveKey=solveKey, evaltmpkw... )[1]

  return getCoordinates(getVariableType(dfg, trgsym), pt)
end


"""
$SIGNATURES

Solve for frontal values only with values in seprarators fixed
  
DevNotes
- WIP
- Relates to: https://github.com/JuliaRobotics/IncrementalInference.jl/issues/466#issuecomment-562556953
- Consolidation
  - Definitely with [`solveFactorParameteric`](@ref)
  - maybe with [`solveGraphParametric`](@ref)
"""
function solveConditionalsParametric(fg::AbstractDFG,
                                    frontals::Vector{Symbol},
                                    separators::Vector{Symbol} = setdiff(listVariables(fg), frontals);
                                    solvekey::Symbol=:parametric,
                                    autodiff = :forward,
                                    algorithm=Optim.BFGS,
                                    algorithmkwargs=(), # add manifold to overwrite computed one
                                    options = Optim.Options(allow_f_increases=true,
                                                            time_limit = 100,
                                                            # show_trace = true,
                                                            # show_every = 1,
                                                            ))

  varIds = [frontals; separators]

  sfg = issetequal(varIds, listVariables(fg)) ? fg : buildSubgraph(fg, varIds, 1)

  flatvar = FlatVariables(fg, varIds)

  for vId in varIds
    p = getVariableSolverData(fg, vId, solvekey).val[1]
    flatvar[vId] =getCoordinates(getVariableType(fg,vId), p)
  end
  initValues = flatvar.X

  frontalsLength = sum(map(v->getDimension(getVariable(fg, v)), frontals))


  # build variables for frontals and seperators
  # fX = view(initValues, 1:frontalsLength)
  fX = initValues[1:frontalsLength]
  # sX = view(initValues, (frontalsLength+1):length(initValues))
  sX = initValues[frontalsLength+1:end]

  alg = algorithm(; algorithmkwargs...)
  # alg = algorithm(; algorithmkwargs...)
  cfd = calcFactorMahalanobisDict(sfg)
  tdtotalCost = Optim.TwiceDifferentiable((x)->_totalCost(fg, cfd, flatvar, [x;sX]), fX, autodiff = autodiff)

  # result = Optim.optimize((x)->_totalCost(fg, flatvar, [x;sX]), fX, alg, options)
  result = Optim.optimize(tdtotalCost, fX, alg, options)

  rv = Optim.minimizer(result)

  H = Optim.hessian!(tdtotalCost, rv)

  Σ = pinv(H)

  d = OrderedDict{Symbol,NamedTuple{(:val, :cov),Tuple{AbstractArray,Matrix{Float64}}}}()

  for key in frontals
    r = flatvar.idx[key]
    p = getPoint(getVariableType(fg, key), rv[r])
    push!(d,key=>(val=p,cov=Σ[r,r]))
  end
  
  return (opti=d, stat=result, varIds=flatvar.idx, Σ=Σ)
end