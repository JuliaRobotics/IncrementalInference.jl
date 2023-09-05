module IncrInfrDiffEqFactorExt

@info "IncrementalInference.jl is loading extensions related to DifferentialEquations.jl"

using DifferentialEquations
import DifferentialEquations: solve

using Dates

using IncrementalInference
import IncrementalInference: getSample, getManifold, DERelative
import IncrementalInference: sampleFactor

using DocStringExtensions

export DERelative



getManifold(de::DERelative{T}) where {T} = getManifold(de.domain)

"""
$SIGNATURES

Calculate a DifferentialEquations.jl ready `tspan::Tuple{Float64,Float64}` from DFGVariables.

DevNotes
- TODO does not yet incorporate Xi.nanosecond field.
- TODO does not handle timezone crossing properly yet.
"""
function _calcTimespan(Xi::AbstractVector{<:DFGVariable})
  #
  tsmps = getTimestamp.(Xi[1:2]) .|> DateTime .|> datetime2unix
  # toffs = (tsmps .- tsmps[1]) .|> x-> elemType(x.value*1e-3)
  return (tsmps...,)
end
# Notes
# - Can change numerical data return type using an additional first argument, `_calcTimespan(Float32, Xi)`.
# _calcTimespan(Xi::AbstractVector{<:DFGVariable}) = _calcTimespan(Float64, Xi)

# performance helper function, FIXME not compatible with all multihypo cases
_maketuplebeyond2args = (w1 = nothing, w2 = nothing, w3_...) -> (w3_...,)

function DERelative(
  Xi::AbstractVector{<:DFGVariable},
  domain::Type{<:InferenceVariable},
  f::Function,
  data = () -> ();
  dt::Real = 1,
  state0::AbstractVector{<:Real} = zeros(getDimension(domain)),
  state1::AbstractVector{<:Real} = zeros(getDimension(domain)),
  tspan::Tuple{<:Real, <:Real} = _calcTimespan(Xi),
  problemType = DiscreteProblem,
)
  #
  datatuple = if 2 < length(Xi)
    datavec = getDimension.([_maketuplebeyond2args(Xi...)...]) .|> x -> zeros(x)
    (data, datavec...)
  else
    data
  end
  # forward time problem
  fproblem = problemType(f, state0, tspan, datatuple; dt = dt)
  # backward time problem
  bproblem = problemType(f, state1, (tspan[2], tspan[1]), datatuple; dt = -dt)
  # build the IIF recognizable object
  return DERelative(domain, fproblem, bproblem, datatuple, getSample)
end

function DERelative(
  dfg::AbstractDFG,
  labels::AbstractVector{Symbol},
  domain::Type{<:InferenceVariable},
  f::Function,
  data = () -> ();
  Xi::AbstractArray{<:DFGVariable} = getVariable.(dfg, labels),
  dt::Real = 1,
  state0::AbstractVector{<:Real} = zeros(getDimension(domain)),
  state1::AbstractVector{<:Real} = zeros(getDimension(domain)),
  tspan::Tuple{<:Real, <:Real} = _calcTimespan(Xi),
  problemType = DiscreteProblem,
)
  return DERelative(
    Xi,
    domain,
    f,
    data;
    dt = dt,
    state0 = state0,
    state1 = state1,
    tspan = tspan,
    problemType = problemType,
  )
end
#
#

# n-ary factor: Xtra splat are variable points (X3::Matrix, X4::Matrix,...)
function _solveFactorODE!(measArr, prob, u0pts, Xtra...)
  # happens when more variables (n-ary) must be included in DE solve
  for (xid, xtra) in enumerate(Xtra)
    # update the data register before ODE solver calls the function
    prob.p[xid + 1][:] = xtra[:]
  end

  # set the initial condition
  prob.u0[:] = u0pts[:]
  sol = DifferentialEquations.solve(prob)

  # extract solution from solved ode
  measArr[:] = sol.u[end]
  return sol
end

# # # output for AbstractRelative is tangents (but currently we working in coordinates for integration with DiffEqs)
# # # FIXME, how to consolidate DERelative with parametric solve which currently only goes through getMeasurementParametric
# function getSample(cf::CalcFactor{<:DERelative})
#   #
#   oder = cf.factor

#   # how many trajectories to propagate?
#   # @show getLabel(cf.fullvariables[2]), getDimension(cf.fullvariables[2])
#   meas = zeros(getDimension(cf.fullvariables[2])) 

#   # pick forward or backward direction
#   # set boundary condition
#   u0pts = if cf.solvefor == 1
#     # backward direction
#     prob = oder.backwardProblem
#     addOp, diffOp, _, _ = AMP.buildHybridManifoldCallbacks(
#       convert(Tuple, getManifold(getVariableType(cf.fullvariables[1]))),
#     )
#     # FIXME use ccw.varValsAll containter?
#     (getBelief(cf.fullvariables[2]) |> getPoints)[cf._sampleIdx]
#   else
#     # forward backward
#     prob = oder.forwardProblem
#     # buffer manifold operations for use during factor evaluation
#     addOp, diffOp, _, _ = AMP.buildHybridManifoldCallbacks(
#       convert(Tuple, getManifold(getVariableType(cf.fullvariables[2]))),
#     )
#     # FIXME use ccw.varValsAll containter?
#     (getBelief(cf.fullvariables[1]) |> getPoints)[cf._sampleIdx]
#   end

#   # solve likely elements
#   # TODO, does this respect hyporecipe ???
#   # TBD check if cf._legacyParams == ccw.varValsAll???
#   idxArr = (k -> cf._legacyParams[k][cf._sampleIdx]).(1:length(cf._legacyParams))
#   _solveFactorODE!(meas, prob, u0pts, _maketuplebeyond2args(idxArr...)...)
#   # _solveFactorODE!(meas, prob, u0pts, i, _maketuplebeyond2args(cf._legacyParams...)...)

#   return meas, diffOp
# end


# NOTE see #1025, CalcFactor should fix `multihypo=` in `cf.__` fields; OBSOLETE
function (cf::CalcFactor{<:DERelative})(measurement, X...)
  #
  meas1 = measurement[1]
  diffOp = measurement[2]

  oderel = cf.factor

  # work on-manifold
  # diffOp = meas[2]
  # if backwardSolve else forward

  # check direction

  solveforIdx = cf.solvefor

  if solveforIdx > 2
    # need to recalculate new ODE (forward) for change in parameters (solving for 3rd or higher variable)
    solveforIdx = 2
    # use forward solve for all solvefor not in [1;2]
    u0pts = getBelief(cf.fullvariables[1]) |> getPoints
    # update parameters for additional variables
    _solveFactorODE!(
      meas1,
      oderel.forwardProblem,
      u0pts[cf._sampleIdx],
      _maketuplebeyond2args(X...)...,
    )
  end

  # find the difference between measured and predicted.
  ## assuming the ODE integrated from current X1 through to predicted X2 (ie `meas1[:,idx]`)
  ## FIXME, obviously this is not going to work for more compilcated groups/manifolds -- must fix this soon!
  # @show cf._sampleIdx, solveforIdx, meas1

  #FIXME 
  res = zeros(size(X[2], 1))
  for i = 1:size(X[2], 1)
    # diffop( test, reference )   <===>   ΔX = test \ reference
    res[i] = diffOp[i](X[solveforIdx][i], meas1[i])
  end
  return res
end




## =========================================================================
## MAYBE legacy

# FIXME see #1025, `multihypo=` will not work properly yet
function IncrementalInference.sampleFactor(cf::CalcFactor{<:DERelative}, N::Int = 1)
  #
  oder = cf.factor

  # how many trajectories to propagate?
  # @show getLabel(cf.fullvariables[2]), getDimension(cf.fullvariables[2])
  meas = [zeros(getDimension(cf.fullvariables[2])) for _ = 1:N]

  # pick forward or backward direction
  # set boundary condition
  u0pts = if cf.solvefor == 1
    # backward direction
    prob = oder.backwardProblem
    addOp, diffOp, _, _ = AMP.buildHybridManifoldCallbacks(
      convert(Tuple, getManifold(getVariableType(cf.fullvariables[1]))),
    )
    getBelief(cf.fullvariables[2]) |> getPoints
  else
    # forward backward
    prob = oder.forwardProblem
    # buffer manifold operations for use during factor evaluation
    addOp, diffOp, _, _ = AMP.buildHybridManifoldCallbacks(
      convert(Tuple, getManifold(getVariableType(cf.fullvariables[2]))),
    )
    getBelief(cf.fullvariables[1]) |> getPoints
  end

  # solve likely elements
  for i = 1:N
    # TODO, does this respect hyporecipe ???
    idxArr = (k -> cf._legacyParams[k][i]).(1:length(cf._legacyParams))
    _solveFactorODE!(meas[i], prob, u0pts[i], _maketuplebeyond2args(idxArr...)...)
    # _solveFactorODE!(meas, prob, u0pts, i, _maketuplebeyond2args(cf._legacyParams...)...)
  end

  return map(x -> (x, diffOp), meas)
end
# getDimension(oderel.domain)





## the function
# ode.problem.f.f

#

end # module