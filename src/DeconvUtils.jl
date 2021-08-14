# series of deconvolution tools


export selectFactorType
export approxDeconv, deconvSolveKey


## Initial version of selecting the dimension of a factor -- will be consolidated with existing infrastructure later


"""
    $SIGNATURES

Inverse solve of predicted noise value and returns tuple of (newly predicted, and known "measured" 
noise) values.

Notes
- Only works for first value in measurement::Tuple at this stage.
- "measured" is used as starting point for the "predicted" values solve.
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
function approxDeconv(fcto::DFGFactor,
                      ccw::CommonConvWrapper = _getCCW(fcto);
                      N::Int=100,
                      measurement::Tuple=sampleFactor(ccw, N),
                      retries::Int=3  )
  #
  # but what if this is a partial factor -- is that important for general cases in deconv?
  _setCCWDecisionDimsConv!(ccw)
  
  # FIXME This does not incorporate multihypo??
  varsyms = getVariableOrder(fcto)
  # vars = getPoints.(getBelief.(dfg, varsyms, solveKey) )
  fcttype = getFactorType(fcto)
  
  # get measurement dimension
  zDim = _getZDim(fcto)
  # TODO consider using ccw.cpt[thrid].res # likely needs resizing
  res_ = zeros(zDim)
  # TODO, consolidate fmd with getSample/sampleFactor and _buildLambda
  fctSmpls = deepcopy(measurement[1])
  fmd = _getFMdThread(ccw)
  
  # TODO assuming vector on only first container in measurement::Tuple
  makeTarget = (i) -> measurement[1][i] # TODO does not support copy-primitive types like Float64, only Ref()
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
  
  for idx in 1:N
    # towards each particle in their own thread (not 100% ready yet, factors should be separate memory)
    thrid = Threads.threadid()
    cpt_ = ccw.cpt[thrid]
    targeti_ = makeTarget(idx)
    
    # TODO must first resolve hypothesis selection before unrolling them -- deferred #1096
    cpt_.activehypo = hyporecipe.activehypo[2][2]

    onehypo!, _ = _buildCalcFactorLambdaSample( ccw,
                                                idx,
                                                cpt_,
                                                targeti_,
                                                measurement  )
    #
    
    # lambda with which to find best measurement values
    hypoObj = (tgt) -> (targeti_.=tgt; onehypo!() )

    # find solution via SubArray view pointing to original memory location
    if fcttype isa AbstractManifoldMinimize
      sfidx = ccw.varidx
      targeti_ .= _solveLambdaNumeric(fcttype, hypoObj, res_, measurement[1][idx], ccw.vartypes[sfidx](), islen1)
    else
      targeti_ .= _solveLambdaNumeric(fcttype, hypoObj, res_, measurement[1][idx], islen1)
    end

  end

  # return (deconv-prediction-result, independent-measurement)
  return measurement[1], fctSmpls
end



"""
    $SIGNATURES

Generalized deconvolution to find the predicted measurement values of the factor `fctsym` in `dfg`.
Inverse solve of predicted noise value and returns tuple of (newly predicted, and known "measured" noise) values.

Notes
- Opposite operation contained in `approxConv`.
- For more notes see [`solveFactorMeasurements`](@ref).

Related

[`approxConv`](@ref), `deconvSolveKey`
"""
function approxDeconv(dfg::AbstractDFG, 
                      fctsym::Symbol, 
                      solveKey::Symbol=:default;
                      retries::Int=3  )
  #
  
  # which factor
  fct = getFactor(dfg, fctsym)
  pts = getPoints(getBelief(dfg,getVariableOrder(fct)[1], solveKey))
  N = length(pts)
  pts = approxDeconv(fct, N=N, retries=retries )
  return pts
end


"""
    $SIGNATURES

Calculate the relative difference between two variables and across solveKeys.

Example

```julia
fg = generateCanonicalFG_Kaess()
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
function deconvSolveKey(dfg::AbstractDFG, 
                        refSym::Symbol, 
                        refKey::Symbol, 
                        tstSym::Symbol, 
                        tstKey::Symbol;  
                        tfg = initfg()  )
  #
  # create a new temporary factor graph for calculations

  # add the first "reference" variable
  Xref = getBelief(dfg, refSym, refKey)
  refSym_ = Symbol(refSym, "_ref")
  refVarType = getVariableType(dfg, refSym)
  if !exists(tfg, refSym_)
    addVariable!(tfg, refSym_, refVarType)
  end
  initManual!(tfg, refSym_, Xref)

  # add the second "test" variable
  tstVarType = getVariableType(dfg, tstSym)
  Xtst = getBelief(dfg, tstSym, tstKey)
  tstSym_ = Symbol(tstSym, "_tst")
  if !exists(tfg, tstSym_)
    addVariable!(tfg, tstSym_, tstVarType)
  end
  initManual!(tfg, tstSym_, Xtst)

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
