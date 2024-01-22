# series of deconvolution tools

## Initial version of selecting the dimension of a factor -- will be consolidated with existing infrastructure later

"""
    $SIGNATURES

Inverse solve of predicted noise value and returns tuple of (newly calculated-predicted, and known measurements) values.

Notes
- Only works for first value in measurement::Tuple at this stage.
- "measured" is used as starting point for the "calculated-predicted" values solve.
- Not all factor evaluation cases are support yet.
- NOTE only works on `.threadid()==1` at present, see #1094
- This function is still part of the initial implementation and needs a lot of generalization improvements.

DevNotes
- TODO Test for various cases with multiple variables.
- TODO make multithread-safe, and able, see #1094
- TODO Test for cases with `nullhypo`
- FIXME FactorMetadata object for all use-cases, not just empty object.
- TODO resolve #1096 (multihypo)
  - TODO Test cases for `multihypo`.
- TODO figure out if there is a way to consolidate with evalFactor and approxConv?
  - basically how to do deconv for just one sample with unique values (wrt TAF)
- TODO N should not be hardcoded to 100

Related

[`approxDeconv`](@ref), [`_solveCCWNumeric!`](@ref)
"""
function approxDeconv(
  fcto::DFGFactor,
  ccw::CommonConvWrapper = _getCCW(fcto);
  N::Int = 100,
  measurement::AbstractVector = sampleFactor(ccw, N),
  retries::Int = 3,
)
  #
  # FIXME needs xDim for all variables at once? xDim = 0 likely to break?

  # but what if this is a partial factor -- is that important for general cases in deconv?
  _setCCWDecisionDimsConv!(ccw, 0) # ccwl.xDim used to hold the last forward solve getDimension(getVariableType(Xi[sfidx]))

  # FIXME This does not incorporate multihypo??
  varsyms = getVariableOrder(fcto)
  # vars = getPoints.(getBelief.(dfg, varsyms, solveKey) )
  fcttype = getFactorType(fcto)

  # get measurement dimension
  zDim = _getZDim(fcto)
  # TODO consider using ccw.cpt[thrid].res # likely needs resizing
  res_ = zeros(zDim)
  # TODO, consolidate fmd with getSample/sampleFactor and _buildLambda
  fctSmpls = deepcopy(measurement)

  # TODO assuming vector on only first container in measurement::Tuple
  makeTarget = (smpidx) -> measurement[smpidx] # TODO does not support copy-primitive types like Float64, only Ref()
  # makeTarget = (i) -> view(measurement[1][i],:)
  # makeTarget = (i) -> view(measurement[1], :, i)

  # NOTE 
  # build a lambda that incorporates the multihypo selections
  # set these first
  # ccw.cpt[].activehypo / .p / .params  # params should already be set from construction
  hyporecipe = _prepareHypoRecipe!(nothing, N, 0, length(varsyms))
  # Juila 1.7 allows destructure assign `(;a,b) = namedtype`
  # certainidx, allelements, activehypo, mhidx = 
  # only doing the current active hypo
  @assert hyporecipe.activehypo[2][1] == 1 "deconv was expecting hypothesis nr == (1, 1:d)"

  islen1 = zDim == 1

  for idx = 1:N
    # towards each particle in their own thread (not 100% ready yet, factors should be separate memory)
    target_smpl = makeTarget(idx)

    # TODO must first resolve hypothesis selection before unrolling them -- deferred #1096
    resize!(ccw.hyporecipe.activehypo, length(hyporecipe.activehypo[2][2]))
    ccw.hyporecipe.activehypo[:] = hyporecipe.activehypo[2][2]

    onehypo! = _buildCalcFactorLambdaSample(ccw, idx, measurement)
    #

    # lambda with which to find best measurement values
    function hypoObj(tgt)
      # copyto!(target_smpl, tgt)
      measurement[idx] = tgt
      return onehypo!()
    end
    # hypoObj = (tgt) -> (target_smpl .= tgt; onehypo!())

    # find solution via SubArray view pointing to original memory location
    if fcttype isa AbstractManifoldMinimize
      error("Fix dispatch on AbstractManifoldMinimize")
    else
      ts = _solveLambdaNumeric(fcttype, hypoObj, res_, measurement[idx], islen1)
      measurement[idx] = ts
    end
  end

  # return (deconv-prediction-result, independent-measurement)
  # r_meas = map(m->m[1], measurement)
  # r_fctSmpls = map(m->m[1], fctSmpls)
  return measurement, fctSmpls
end

# TBD deprecate use of xDim
function approxDeconv(
  fcto::DFGFactor{<:CommonConvWrapper{<:AbstractManifoldMinimize}},
  ccw::CommonConvWrapper = _getCCW(fcto);
  N::Int = 100,
  measurement::AbstractVector = sampleFactor(ccw, N),
  retries=nothing,
)
  if !isnothing(retries)
    Base.depwarn(
      "approxDeconv kwarg retries is not used",
      :approxDeconv,
    )
  end
  # but what if this is a partial factor -- is that important for general cases in deconv?
  _setCCWDecisionDimsConv!(ccw, 0)

  varsyms = getVariableOrder(fcto)
  
  # TODO assuming vector on only first container in measurement::Tuple  # TBD How should user dispatch fancy tuple measurements on deconv.
  
  # NOTE 
  # build a lambda that incorporates the multihypo selections
  # deconv has to solve for the best matching for particles
  # FIXME This does not incorporate multihypo, Apply hyporecipe to full variable order list. But remember hyporecipe assignment must be found (NPhard)
  hyporecipe = _prepareHypoRecipe!(nothing, N, 0, length(varsyms))
  # only doing the current active hypo
  @assert hyporecipe.activehypo[2][1] == 1 "deconv was expecting hypothesis nr == (1, 1:d)"

  # get measurement dimension
  zDim = _getZDim(fcto)
  islen1 = zDim == 1

  #make a copy of the original measurement before mutating it
  sampled_meas = deepcopy(measurement)

  fcttype = getFactorType(fcto)
  
  for idx = 1:N

    # TODO must first resolve hypothesis selection before unrolling them -- deferred #1096
    resize!(ccw.hyporecipe.activehypo, length(hyporecipe.activehypo[2][2]))
    ccw.hyporecipe.activehypo[:] = hyporecipe.activehypo[2][2]
    #TODO why is this resize in the loop?

    # Create a CalcFactor functor of the correct hypo.
    _hypoCalcFactor = _buildHypoCalcFactor(ccw, idx)

    ts = _solveLambdaNumericMeas(fcttype, _hypoCalcFactor, measurement[idx], islen1)
    measurement[idx] = ts

  end

  return measurement, sampled_meas
end

"""
    $SIGNATURES

Generalized deconvolution to find the predicted measurement values of the factor `fctsym` in `dfg`.
Inverse solve of predicted noise value and returns tuple of (newly predicted, and known "measured" noise) values.

Notes
- Opposite operation contained in `approxConvBelief`.
- For more notes see [`solveFactorMeasurements`](@ref).

Related

[`approxConvBelief`](@ref), `deconvSolveKey`
"""
function approxDeconv(
  dfg::AbstractDFG,
  fctsym::Symbol,
  solveKey::Symbol = :default;
  retries::Int = 3,
)
  #

  # which factor
  fct = getFactor(dfg, fctsym)
  pts = getPoints(getBelief(dfg, getVariableOrder(fct)[1], solveKey))
  N = length(pts)
  pts = approxDeconv(fct; N = N, retries = retries)
  return pts
end

function approxDeconv(
  dfg::AbstractDFG,
  fctlbl::Symbol,
  factorType::AbstractRelative,
  solveKey::Symbol = :default;
  tfg::AbstractDFG = initfg(),
  retries::Int = 3,
)
  #

  # build a local temporary graph copy containing the same values but user requested factor type.
  fct = getFactor(dfg, fctlbl)
  fctT = getFactorType(fct)
  lbls = getVariableOrder(fct)
  for lb in lbls
    exists(tfg, lb) ? nothing : addVariable!(tfg, lb, getVariableType(dfg, lb))
    initVariable!(tfg, lb, getBelief(dfg, lb, solveKey))
  end

  # add factor type requested by user
  f_ = addFactor!(tfg, lbls, factorType; graphinit = false)

  # peform the deconvolution operation on the temporary graph with user desired factor instead.
  return approxDeconv(tfg, getLabel(f_); retries = retries)
end

# try default constructor
function approxDeconv(
  dfg::AbstractDFG,
  fctlbl::Symbol,
  factorType::Type{<:AbstractRelative},
  w...;
  kw...,
)
  return approxDeconv(dfg, fctlbl, factorType(), w...; kw...)
end
#

function approxDeconvBelief(dfg::AbstractDFG, lb::Symbol, w...; kw...)
  return manikde!(
    getManifold(getFactorType(dfg, lb)),
    approxDeconv(dfg, lb, w...; kw...)[1],
  )
end

"""
    $SIGNATURES

Calculate the relative difference between two variables and across solveKeys.

Example

```julia
fg = generateGraph_Kaess()
solveTree!(fg, storeOld=true)
# calculate the relative motion induced by the solver from init to solve.
pts = deconvSolveKey(fg, :x1, :default, :x1, :graphinit)
```

Notes
- Can pass user `tfg::AbstractDFG` for better in-place operations.

DevNotes
- TODO use dfg, rather than building new tfg internally.

Related

[`approxDeconv`](@ref), [`mmd`](@ref)
"""
function deconvSolveKey(
  dfg::AbstractDFG,
  refSym::Symbol,
  refKey::Symbol,
  tstSym::Symbol,
  tstKey::Symbol;
  tfg = initfg(),
)
  #
  # create a new temporary factor graph for calculations

  # add the first "reference" variable
  Xref = getBelief(dfg, refSym, refKey)
  refSym_ = Symbol(refSym, "_ref")
  refVarType = getVariableType(dfg, refSym)
  if !exists(tfg, refSym_)
    addVariable!(tfg, refSym_, refVarType)
  end
  initVariable!(tfg, refSym_, Xref)

  # add the second "test" variable
  tstVarType = getVariableType(dfg, tstSym)
  Xtst = getBelief(dfg, tstSym, tstKey)
  tstSym_ = Symbol(tstSym, "_tst")
  if !exists(tfg, tstSym_)
    addVariable!(tfg, tstSym_, tstVarType)
  end
  initVariable!(tfg, tstSym_, Xtst)

  # add the new dummy factor with default manifold for computations
  fctType = selectFactorType(refVarType, tstVarType)
  nf = addFactor!(tfg, [refSym_; tstSym_], fctType())

  # TODO connect from dfg all other data that might form part of FactorMetadata in tfg
  pts = approxDeconv(tfg, nf.label)

  # assuming tfg was passed in by the user
  deleteFactor!(tfg, nf.label)

  # return result
  return pts, fctType
end

#
