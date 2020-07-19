# series of deconvolution tools


export selectFactorType
export buildFactorDefault
export solveFactorMeasurements
export buildGraphLikelihoodsDifferential!


## Initial version of selecting the dimension of a factor -- will be consolidated with existing infrastructure later

getManifolds(vartype::Type{LinearConditional}) = (:Euclid,)


##

"""
    $SIGNATURES

First hacky version to return which factor type to use between two variables of types T1 and T2.
"""
selectFactorType(T1::Type{ContinuousScalar}, T2::Type{ContinuousScalar}) = LinearConditional
selectFactorType(T1::InferenceVariable, T2::InferenceVariable) = selectFactorType(typeof(T1), typeof(T2))
selectFactorType(dfg::AbstractDFG, s1::Symbol, s2::Symbol) = selectFactorType( getVariableType(dfg, s1), getVariableType(dfg, s2) )


"""
    $SIGNATURES

Inverse solve of predicted noise value and returns tuple of (newly predicted, and known "measured" noise) values.

Notes
- "measured" is used as starting point for the "predicted" values solve.

DevNotes
- This function is still part of the initial implementation and needs a lot of generalization improvements.
"""
function solveFactorMeasurements(dfg::AbstractDFG,
                                 fctsym::Symbol  )
  #
  fcto = getFactor(dfg, fctsym)
  varsyms = getVariableOrder(fcto)
  vars = map(x->getPoints(getKDE(dfg,x)), varsyms)
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
      if isa(fcttype, FunctorPairwiseMinimize)
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
      elseif isa(fcttype, FunctorPairwise)
        ggnl = (rs, dm) -> fcttype(rs,ud,idx,makemeas!(idx, meas, dm),vars...)
        r = nlsolve(ggnl, meas[1][:,idx])
        break
      elseif isa(fcttype, FunctorSingleton)
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




function approxDeconv(dfg::AbstractDFG, sym1::Symbol, sym2::Symbol)
  # which factor
  fnames = intersect(ls(dfg, sym1), ls(dfg,sym2))

  @assert length(fnames) == 1 "approxDeconv cannot yet handle multiple parallel factors between the same variables -- this is a TODO, please open an issue with IncrementalInference"
  pts = solveFactorMeasurements(dfg, fnames[1])


end


"""
    $SIGNATURES

Build from a `LikelihoodMessage` a temporary distributed factor graph object containing differential
information likelihood factors based on values in the messages.

Notes
- Modifies tfg argument by adding `:UPWARD_DIFFERENTIAL` factors.

DevNotes
- Initial version which only works for Pose2 and Point2 at this stage.
"""
function buildGraphLikelihoodsDifferential!(msgs::LikelihoodMessage,
                                            tfg=initfg() )
  # create new local dfg and add all the variables with data
  for (label, val) in msgs.belief
    addVariable!(tfg, label, val.softtype)
    initManual!(tfg, label, manikde!(val))
  end

  # list all variables in order of dimension size
  alreadylist = Symbol[]
  listVarByDim = ls(tfg)
  listDims = getDimension.(getVariable.(tfg,listVarByDim))
  per = sortperm(listDims, rev=true)
  listVarDec = listVarByDim[per]
  listVarAcc = reverse(listVarDec)
  # add all differential factors (without deconvolution values)
  for sym1_ in listVarDec
    push!(alreadylist, sym1_)
    for sym2_ in setdiff(listVarAcc, alreadylist)
      nfactype = selectFactorType(tfg, sym1_, sym2_)
      # assume default helper function # buildFactorDefault(nfactype)
      nfct = nfactype()
      afc = addFactor!(tfg, [sym1_;sym2_], nfct, graphinit=false, tags=[:DUMMY;])
      # calculate the general deconvolution between variables
      pts = solveFactorMeasurements(tfg, afc.label)
      newBel = manikde!(pts[1], getManifolds(nfactype))
      # replace dummy factor with real deconv factor using manikde approx belief measurement
      fullFct = nfactype(newBel)
      deleteFactor!(tfg, afc.label)
      addFactor!( tfg, [sym1_;sym2_], fullFct, graphinit=false, tags=[:LIKELIHOODMESSAGE; :UPWARD_DIFFERENTIAL] )
    end
  end

  return tfg
end


#
