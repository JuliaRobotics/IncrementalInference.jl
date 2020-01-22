# tether utils

export cont2disc
export rebaseFactorVariable!
export getFactorMean
export solveBinaryFactorParameteric, accumulateFactorMeans

"""
    $SIGNATURES

Standard mean and covariance propagation for linear systems.  Directly equivalent to Kalman filtering methods.

Notes
- Does the proper continuous (Qc) to discrete process noise (Qd) calculation -- as per Farrell, 2008.
- Used downstream for in-time Gaussian mean and covariance propagation.
"""
function cont2disc(F::Matrix{Float64},
                   G::Matrix{Float64},
                   Qc::Matrix{Float64},
                   dt::Float64,
                   Phik::Matrix{Float64}=Matrix{Float64}(LinearAlgebra.I, 0,0) )
    #
    fr,fc = size(F)
    gr,gc = size(G)

    # two bits of new memory allocated
    M1 = zeros(fc+gc,fc+gc)
    M2 = zeros(fr+fc,fr+fc)

    M1[1:fr,1:fc] = F
    M1[1:gr,(fc+1):end] = G #1:gr,(fc+1):(fc+gc)

    # must convert to propagateLinSystem call, use trapezoidal
    Md1 = exp(M1*dt) # heavy lifting here
    Phi = size(Phik,1) == 0 ? Md1[1:fr,1:fc] : Phik
    Gamma = Md1[1:fr,(fc+1):end]

    #M2 = [[-F';(G*Qc*G')']';[zeros(9,9);F]'] # easy concat
    GQG = (G*Qc*G')
    gqgr, gqgc = size(GQG)
    M2[1:fc,1:fr] = -F
    M2[1:fr,(fc+1):end] = GQG
    M2[(fr+1):end,(fc+1):end] = F'

    Md2 = exp(M2*dt) # heavy lifting here
    Qd = Phi * Md2[1:fr,(fc+1):end] #Qd = Phi*(Md2[1:fr,(fc+1):end])

    # Qd = GQG*dt;

    return Phi, Gamma, Qd
end



"""
    $SIGNATURES

Helper function to modify factor connectivity to variables.

Notes
- Developed for updating a dead reckoning odometry factor.
- Arguments are order sensitive.
"""
function rebaseFactorVariable!(dfg::AbstractDFG,
                               fctsym::Symbol,
                               newvars::Vector{Symbol};
                               rmDisconnected::Bool=true,
                               autoinit::Bool=false  )::Nothing
  #
  # check that all new variables are available
  @assert sum(map(x->exists(dfg, x), newvars)) == length(newvars)

  # get existing factor details
  fct = getFactor(dfg, fctsym)
  fcttype = getFactorType(fct)
  mh = getMultihypoDistribution(fct)

  # get old vars
  oldvars = getVariableOrder(fct)

  # delete old factor from graph
  deleteFactor!(dfg, fctsym)

  # add the factor back into graph against new variables
  addFactor!(dfg, newvars, fcttype, autoinit=autoinit, multihypo=mh)

  # clean up disconnected variables if requested
  if rmDisconnected
    for ov in oldvars
      # find variables that are not connected to anything
      if length(ls(dfg, ov)) == 0
        deleteVariable!(dfg, ov)
      end
    end
  end

  return nothing
end


"""
    $SIGNATURES

Recover the mean (Gaussian) or estimate stochastic mean (non-Gaussian) value stored in a factor measurement.

Related

accumulateFactorMeans, solveBinaryFactorParameteric
"""
function getFactorMean(fct::FunctorInferenceType)
  fctt = typeof(getFactorType(fct))
  error("no getFactorMean defined for $(fctt.name), has fields $(fieldnames(fctt))")
end

getFactorMean(fct::Normal) = fct.μ
getFactorMean(fct::MvNormal) = fct.μ
getFactorMean(fct::BallTreeDensity) = getKDEMean(fct)
getFactorMean(fct::AliasingScalarSampler) = Statistics.mean(rand(fct,1000))

getFactorMean(fct::DFGFactor) = getFactorMean(getFactorType(fct))

getFactorMean(dfg::AbstractDFG, fctsym::Symbol) = getFactorMean(getFactor(dfg, fctsym))

"""
    $SIGNATURES

Helper function to propagate a parametric estimate along a factor chain.

Notes
- Not used during mmisam inference.
- Expected uses are for user analysis of factors and estimates.
- real-time dead reckoning chain prediction.

Related:

approxConv, accumulateFactorMeans, MutablePose2Pose2Gaussian
"""
function solveBinaryFactorParameteric(dfg::AbstractDFG,
                                      fct::DFGFactor,
                                      currval::Vector{Float64},
                                      srcsym::Symbol,
                                      trgsym::Symbol  )::Vector{Float64}
  #
  outdims = getVariableDim(getVariable(dfg, trgsym))
  meas = getFactorType(fct)
  mea = getFactorMean(fct)
  measT = (reshape(mea,:,1),)

  # calculate the projection
  varmask = (1:2)[getVariableOrder(fct) .== trgsym][1]
  pts = approxConvBinary( reshape(currval,:,1), meas, outdims, measT, varidx=varmask )

  # return the result
  @assert length(pts) == outdims
  return pts[:]
end

"""
    $SIGNATURES

Accumulate chains of binary factors---potentially starting from a prior---as a parameteric mean value only.

Notes
- Not used during mmisam inference.
- Expected uses are for user analysis of factors and estimates.
- real-time dead reckoning chain prediction.

Related:

approxConv, solveBinaryFactorParameteric, MutablePose2Pose2Gaussian
"""
function accumulateFactorMeans(dfg::AbstractDFG, fctsyms::Vector{Symbol})

  ## get the starting estimate
  val = zeros(0)
  nextidx = 1
  onePrior = false
  currsym = :null
  if isPrior(dfg, fctsyms[nextidx])
    # if first factor is prior
    @assert !onePrior
    onePrior = true
    val = getFactorMean(dfg, fctsyms[nextidx])
    currsym = ls(dfg, fctsyms[nextidx])[1]
    nextidx += 1
  else
    # get first value from current variable estimate
    vars = getVariableOrder(dfg, fctsyms[nextidx])
    nextsym = 1 < length(fctsyms) ? intersect( vars, ls(dfg, fctsyms[nextidx+1]) ) : vars[end]
    currsym = 1 < length(fctsyms) ? setdiff(vars, nextsym)[1] : vars[1]
    val = calcVariablePPE(dfg, currsym).suggested
  end

  srcsym = currsym
  # Propagate the parametric value along the factor chain
  for fct in map(x->getFactor(dfg, x), fctsyms[nextidx:end])
    # first find direction of solve
    vars = getVariableOrder(fct)
    trgsym = setdiff(vars, [srcsym])[1]
    # varmask = (1:2)[getVariableOrder(fct) .== trgsym][1]
    val = solveBinaryFactorParameteric(dfg,fct,val,srcsym,trgsym)
    srcsym = trgsym
  end

  return val
end
