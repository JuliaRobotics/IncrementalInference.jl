# series of deconvolution tools

export solveFactorMeasurements



"""
    $SIGNATURES

Inverse solve of predicted noise value and returns the associated "measured" noise value (also used as starting point for the solve).
"""
function solveFactorMeasurements(dfg::AbstractDFG,
                                 fctsym::Symbol  )
  #
  fcto = getFactor(dfg, fctsym)
  varsyms = fcto._variableOrderSymbols
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

  var1 = getVariable(dfg, sym1)
  var2 = getVariable(dfg, sym2)




end




#
