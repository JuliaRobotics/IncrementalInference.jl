# series of deconvolution tools


export selectFactorType
export solveFactorMeasurements
export approxDeconv, deconvSolveKey


## Initial version of selecting the dimension of a factor -- will be consolidated with existing infrastructure later

getDomain(::InstanceType{LinearConditional}) = ContinuousScalar

getManifolds(fctType::Type{LinearConditional}) = getManifolds(getDomain(fctType))


##

"""
    $SIGNATURES

First hacky version to return which factor type to use between two variables of types T1 and T2.
"""
selectFactorType(T1::Type{ContinuousScalar}, T2::Type{ContinuousScalar}) = LinearConditional
selectFactorType(T1::Type{<:ContinuousEuclid{N}}, T2::Type{<:ContinuousEuclid{N}}) where N = LinearConditional{N}
selectFactorType(T1::InferenceVariable, T2::InferenceVariable) = selectFactorType(typeof(T1), typeof(T2))
selectFactorType(dfg::AbstractDFG, s1::Symbol, s2::Symbol) = selectFactorType( getVariableType(dfg, s1), getVariableType(dfg, s2) )


"""
    $SIGNATURES

Inverse solve of predicted noise value and returns tuple of (newly predicted, and known "measured" noise) values.

Notes
- "measured" is used as starting point for the "predicted" values solve.
- Not all factor evaluation cases are support yet.

DevNotes
- This function is still part of the initial implementation and needs a lot of generalization improvements.
- FIXME FactorMetadata object for all use-cases, not just empty object.
- Test for various cases with multiple variables.
- Test for cases with `nullhypo` and `multihypo`.
"""
function solveFactorMeasurements(dfg::AbstractDFG,
                                 fctsym::Symbol,
                                 solveKey::Symbol=:default  )
  #
  fcto = getFactor(dfg, fctsym)
  varsyms = getVariableOrder(fcto)
  vars = map(x->getPoints(getBelief(dfg,x,solveKey)), varsyms)
  fcttype = getFactorType(fcto)

  N = size(vars[1])[2]
  ud = FactorMetadata()
  meas = getSample(fcttype, N)
  meas0 = deepcopy(meas[1])
  # get measurement dimension
  zDim = getSolverData(fcto).fnc.zDim
  res = zeros(zDim)

  function makemeas!(i, meas, dm)
    meas[1][:,i] = dm
    return meas
  end

  ggo = (i, dm) -> fcttype(res,ud,i,makemeas!(i, meas, dm),vars...)
  # ggo(1, [0.0;0.0])

  for idx in 1:N
    retry = 10
    while 0 < retry
      if isa(fcttype, AbstractRelativeFactorMinimize)
        r = optimize((x) -> ggo(idx,x), meas[1][:,idx]) # zeros(zDim)
        retry -= 1
        if !r.g_converged
          nsm = getSample(fcttype, 1)
          for count in 1:length(meas)
            meas[count][:,idx] = nsm[count][:,idx]
          end
        else
          break
        end
      elseif isa(fcttype, AbstractRelativeFactor)
        ggnl = (rs, dm) -> fcttype(rs,ud,idx,makemeas!(idx, meas, dm),vars...)
        r = nlsolve(ggnl, meas[1][:,idx])
        break
      elseif isa(fcttype, AbstractPrior)
        # assuming no partials at this point
        meas[1][:,:] .= vars[1][:,:]
        break
      end
    end
    # @assert meas[1][:,idx] == r.minimizer
  end

  # Gadfly.plot(z=(x,y)->ggo(1,[x;y]), xmin=[-pi],xmax=[pi],ymin=[-100.0],ymax=[100.0], Geom.contour)
  return meas[1], meas0
end



"""
    $SIGNATURES

Generalized deconvolution to find the predicted measurement values of the factor `fctsym` in `dfg`.

Notes
- Opposite operation contained in `approxConv`

Related

solveFactorMeasurements, deconvSolveKey, approxConv
"""
function approxDeconv(dfg::AbstractDFG, fctsym::Symbol, solveKey::Symbol=:default)
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

Related

approxDeconv, mmd
"""
function deconvSolveKey(dfg::AbstractDFG, 
                        refSym::Symbol, 
                        refKey::Symbol, 
                        tstSym::Symbol, 
                        tstKey::Symbol  )
  #
  # create a new temporary factor graph for calculations
  tfg = initfg()

  # add the first "reference" variable
  refVarType = getVariableType(dfg, refSym)
  Xref = getBelief(dfg, refSym, refKey)
  refSym_ = Symbol(refSym, "_ref")
  addVariable!(tfg, refSym_, refVarType)
  initManual!(tfg, refSym_, Xref)

  # add the second "test" variable
  tstVarType = getVariableType(dfg, tstSym)
  Xtst = getBelief(dfg, tstSym, tstKey)
  tstSym_ = Symbol(tstSym, "_tst")
  addVariable!(tfg, tstSym_, tstVarType)
  initManual!(tfg, tstSym_, Xtst)

  # add the new dummy factor with default manifold for computations
  fctType = selectFactorType(refVarType, tstVarType)
  nf = addFactor!(tfg, [refSym_; tstSym_], fctType())

  # TODO connect from dfg all other data that might form part of FactorMetadata in tfg
  pts = approxDeconv(tfg, nf.label)

  # return result
  return pts, fctType
end


function mmd(p1::AbstractMatrix{<:Real}, 
             p2::AbstractMatrix{<:Real}, 
             varType::Union{InstanceType{InferenceVariable},InstanceType{FunctorInferenceType}};
             bw::AbstractVector{<:Real}=[0.001;] )
  #
  manis = convert(AMP.Manifold, varType)
  mmd(p1, p2, manis, bw=bw)  
end

function mmd(p1::BallTreeDensity, 
             p2::BallTreeDensity, 
             nodeType::Union{InstanceType{InferenceVariable},InstanceType{FunctorInferenceType}};
             bw::AbstractVector{<:Real}=[0.001;])
  #
  mmd(getPoints(p1), getPoints(p2), nodeType, bw=bw)
end



#
