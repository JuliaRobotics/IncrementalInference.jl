# series of deconvolution tools


export selectFactorType
export solveFactorMeasurements
export approxDeconv, deconvSolveKey


## Initial version of selecting the dimension of a factor -- will be consolidated with existing infrastructure later


"""
    $SIGNATURES

Inverse solve of predicted noise value and returns tuple of (newly predicted, and known "measured" noise) values.

Notes
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

Related

[`approxDeconv`](@ref)
"""
function solveFactorMeasurements( dfg::AbstractDFG,
                                  fctsym::Symbol,
                                  solveKey::Symbol=:default  )
  #
  thrid = Threads.threadid()
  fcto = getFactor(dfg, fctsym)
  ccw = _getCCW(fcto)
  fmd = _getFMdThread(ccw)

  # FIXME This does not incorporate multihypo??
  varsyms = getVariableOrder(fcto)
  vars = map(x->getPoints(getBelief(dfg,x,solveKey)), varsyms)
  fcttype = getFactorType(fcto)

  N = size(vars[1])[2]
  
  # TODO, consolidate this fmd with getSample/freshSamples and _buildLambda

  meas = freshSamples(fcttype, N, fmd) # getSample(fcttype, N)
  meas0 = deepcopy(meas[1])
  # get measurement dimension
  zDim = _getZDim(fcto)

  # TODO consider using ccw.cpt[thrid].res # likely needs resizing
  res_ = zeros(zDim)

  # TODO assuming vector on only first container in meas::Tuple
  makeTarget = (i) -> view(meas[1], :, i)
  
  for idx in 1:N
    targeti_ = makeTarget(idx)

    # NOTE Deferred as #1096
    # TODO must first resolve hypothesis selection before unrolling them
    # build a lambda that incorporates the multihypo selections
    # unrollHypo, _ = _buildCalcFactorLambdaSample( ccw,
    #                                               idx,
    #                                               ccw.cpt[thrid],
    #                                               targeti_,
    #                                               meas  )
    #
    
    # ggo = (res, dm) -> (targeti_.=dm; unrollHypo(res))
    ggo = (res, dm) -> (targeti_.=dm; fcttype(res, fmd, idx, meas, vars...))

    retry = 10
    # FIXME remove the retry steps
    while 0 < retry
      if isa(fcttype, AbstractRelativeMinimize)
        r = if size(targeti_,1) == 1
          optimize((x) -> ggo(res_, x), meas[1][:,idx], BFGS() )
        else
          optimize((x) -> ggo(res_, x), meas[1][:,idx])
        end
        retry -= 1
        if !r.g_converged
          nsm = getSample(fcttype, 1)
          for count in 1:length(meas)
            meas[count][:,idx] = nsm[count][:,1]
          end
        else
          break
        end
      elseif isa(fcttype, AbstractRelativeRoots)
        r = nlsolve(ggo, meas[1][:,idx])
        break
      elseif isa(fcttype, AbstractPrior)
        # return trivial case of meas == meas0
        break
      end
    end
  end

  # return (deconv-prediction-result, independent-measurement)
  return meas[1], meas0
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
                      solveKey::Symbol=:default)
  # which factor
  pts = solveFactorMeasurements(dfg, fctsym, solveKey)
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


"""
    $TYPEDSIGNATURES

Calculate the Kernel Embedding MMD 'distance' between sample points (or kernel density estimates).

Notes
- `bw::Vector=[0.001;]` controls the mmd kernel bandwidths.

Related

`KDE.kld`
"""
function mmd( p1::AbstractMatrix{<:Real}, 
              p2::AbstractMatrix{<:Real}, 
              varType::Union{InstanceType{InferenceVariable},InstanceType{FunctorInferenceType}};
              bw::AbstractVector{<:Real}=[0.001;] )
  #
  manis = convert(AMP.Manifold, varType)
  mmd(p1, p2, manis, bw=bw)  
end

function mmd( p1::BallTreeDensity, 
              p2::BallTreeDensity, 
              nodeType::Union{InstanceType{InferenceVariable},InstanceType{FunctorInferenceType}};
              bw::AbstractVector{<:Real}=[0.001;])
  #
  mmd(getPoints(p1), getPoints(p2), nodeType, bw=bw)
end



#
