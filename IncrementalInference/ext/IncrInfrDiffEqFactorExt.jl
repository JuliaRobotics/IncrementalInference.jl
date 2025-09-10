module IncrInfrDiffEqFactorExt

@info "IncrementalInference.jl is loading extensions related to DifferentialEquations.jl"

import Base: show

using DifferentialEquations
import DifferentialEquations: solve

using Dates

using IncrementalInference
import IncrementalInference: DERelative, _solveFactorODE!
import IncrementalInference: getSample, sampleFactor, getManifold

using DocStringExtensions

export DERelative

using LieGroups: allocate, compose, hat, Identity, vee, log, LieAlgebra


getManifold(de::DERelative{T}) where {T} = getManifold(de.domain)


function Base.show(
  io::IO, 
  ::Union{<:DERelative{T,O},Type{<:DERelative{T,O}}}
) where {T,O}
  println(io, "  DERelative{")
  println(io, "    ", T)
  println(io, "    ", O.name.name)
  println(io, "  }")
  nothing
end

Base.show(
  io::IO, 
  ::MIME"text/plain", 
  der::DERelative
) = show(io, der)


"""
$SIGNATURES

Calculate a DifferentialEquations.jl ready `tspan::Tuple{Float64,Float64}` from DFGVariables.

DevNotes
- TODO does not yet incorporate Xi.nanosecond field.
- TODO does not handle timezone crossing properly yet.
"""
function _calcTimespan(
  Xi::AbstractVector{<:VariableCompute}
)
  #
  tsmps = getTimestamp.(Xi[1:2]) .|> DateTime .|> datetime2unix
  # toffs = (tsmps .- tsmps[1]) .|> x-> elemType(x.value*1e-3)
  return (tsmps...,)
end
# Notes
# - Can change numerical data return type using an additional first argument, `_calcTimespan(Float32, Xi)`.
# _calcTimespan(Xi::AbstractVector{<:VariableCompute}) = _calcTimespan(Float64, Xi)

# performance helper function, FIXME not compatible with all multihypo cases
_maketuplebeyond2args = (w1 = nothing, w2 = nothing, w3_...) -> (w3_...,)

function DERelative(
  Xi::AbstractVector{<:VariableCompute},
  domain::Type{<:StateType},
  f::Function,
  data = () -> ();
  dt::Real = 1,
  state0::AbstractVector{<:Real} = allocate(getPointIdentity(domain)), # zeros(getDimension(domain)),
  state1::AbstractVector{<:Real} = allocate(getPointIdentity(domain)), # zeros(getDimension(domain)),
  tspan::Tuple{<:Real, <:Real} = _calcTimespan(Xi),
  problemType = ODEProblem, # DiscreteProblem,
  keepSolution::Bool = false
)
  #
  datatuple = if 2 < length(Xi)
    datavec = getDimension.([_maketuplebeyond2args(Xi...)...]) .|> x -> zeros(x)
    (data, datavec...)
  else
    data
  end
  # forward time problem
  fproblem = problemType(f, state0, tspan, datatuple; dt)
  # backward time problem
  bproblem = problemType(f, state1, (tspan[2], tspan[1]), datatuple; dt = -dt)

  _keepSolution = if keepSolution
    Vector{typeof((;t=Vector{Float64}(),u=Vector{Vector{Float64}}()))}()
  else
    nothing
  end 

  # build the IIF recognizable object
  return DERelative(
    domain, 
    fproblem, 
    bproblem, 
    datatuple, 
    _keepSolution
  )
end

function DERelative(
  dfg::AbstractDFG,
  labels::AbstractVector{Symbol},
  domain::Type{<:StateType},
  f::Function,
  data = () -> ();
  Xi::AbstractArray{<:VariableCompute} = getVariable.(dfg, labels),
  dt::Real = 1,
  state1::AbstractVector{<:Real} = allocate(getPointIdentity(domain)), #zeros(getDimension(domain)),
  state0::AbstractVector{<:Real} = allocate(getPointIdentity(domain)), #zeros(getDimension(domain)),
  tspan::Tuple{<:Real, <:Real} = _calcTimespan(Xi),
  problemType = DiscreteProblem,
  keepSolution::Bool = false
)
  return DERelative(
    Xi,
    domain,
    f,
    data;
    dt,
    state0,
    state1,
    tspan,
    problemType,
    keepSolution,
  )
end
#
#

# n-ary factor: Xtra splat are variable points (X3::Matrix, X4::Matrix,...)
function _solveFactorODE!(
  measArr, 
  prob, 
  u0pts, 
  Xtra...
)
  # happens when more variables (n-ary) must be included in DE solve
  for (xid, xtra) in enumerate(Xtra)
    # update the data register before ODE solver calls the function
    prob.p[xid + 1][:] = xtra[:] # FIXME, unlikely to work with ArrayPartition, maybe use MArray and `.=`
  end

  # set the initial condition
  prob.u0 .= u0pts

  sol = DifferentialEquations.solve(prob)

  # extract solution from solved ode
  measArr[:] = sol.u[end]
  return sol
end

# # # output for AbstractRelativeObservation is tangents (but currently we working in coordinates for integration with DiffEqs)
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
function (cf::CalcFactor{<:DERelative})(
  measurement, 
  X...
)
  #
  # numerical measurement values
  meas1 = measurement[1]
  # work on-manifold via sampleFactor piggy back of particular manifold definition
  M = measurement[2]
  # lazy factor pointer
  oderel = cf.factor
  # check direction
  solveforIdx = cf.solvefor
  
  # if backwardSolve else forward
  if solveforIdx > 2
    # need to recalculate new ODE (forward) for change in parameters (solving for 3rd or higher variable)
    solveforIdx = 2
    # use forward solve for all solvefor not in [1;2]
    # u0pts = getBelief(cf.fullvariables[1]) |> getPoints
    # update parameters for additional variables
    _solveFactorODE!(
      meas1,
      oderel.forwardProblem,
      X[1], # u0pts[cf._sampleIdx],
      _maketuplebeyond2args(X...)...,
    )
  end

  # find the difference between measured and predicted.
  # assuming the ODE integrated from current X1 through to predicted X2 (ie `meas1[:,idx]`)
  res_ = compose(M, inv(M, X[solveforIdx]), meas1)
  res = vee(LieAlgebra(M), log(M, res_))

  return res
end


# # FIXME see #1025, `multihypo=` will not work properly yet
# function getSample(cf::CalcFactor{<:DERelative})

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
#     cf._legacyParams[2]
#   else
#     # forward backward
#     prob = oder.forwardProblem
#     # buffer manifold operations for use during factor evaluation
#     addOp, diffOp, _, _ = AMP.buildHybridManifoldCallbacks(
#       convert(Tuple, getManifold(getVariableType(cf.fullvariables[2]))),
#     )
#     cf._legacyParams[1]
#   end

#   i = cf._sampleIdx
#   # solve likely elements
#   # TODO, does this respect hyporecipe ???
#   idxArr = (k -> cf._legacyParams[k][i]).(1:length(cf._legacyParams))
#   _solveFactorODE!(meas, prob, u0pts[i], _maketuplebeyond2args(idxArr...)...)
#   # _solveFactorODE!(meas, prob, u0pts, i, _maketuplebeyond2args(cf._legacyParams...)...)

#   return meas, diffOp
# end




## =========================================================================
## MAYBE legacy

# FIXME see #1025, `multihypo=` will not work properly yet
function IncrementalInference.sampleFactor(cf::CalcFactor{<:DERelative}, N::Int = 1)
  #
  oder = cf.factor

  # how many trajectories to propagate?
  # 
  v2T = getVariableType(cf.fullvariables[2])
  meas = [allocate(getPointIdentity(v2T)) for _ = 1:N]
  # meas = [zeros(getDimension(cf.fullvariables[2])) for _ = 1:N]

  # pick forward or backward direction
  # set boundary condition
  u0pts, M = if cf.solvefor == 1
    # backward direction
    prob = oder.backwardProblem
    M_ = getManifold(getVariableType(cf.fullvariables[1]))
    addOp, diffOp, _, _ = AMP.buildHybridManifoldCallbacks(
      AMP._manifoldtuple(M_),
    )
    # getBelief(cf.fullvariables[2]) |> getPoints
    cf._legacyParams[2], M_
  else
    # forward backward
    prob = oder.forwardProblem
    M_ = getManifold(getVariableType(cf.fullvariables[2]))
    # buffer manifold operations for use during factor evaluation
    addOp, diffOp, _, _ = AMP.buildHybridManifoldCallbacks(
      AMP._manifoldtuple(M_),
    )
    # getBelief(cf.fullvariables[1]) |> getPoints
    cf._legacyParams[1], M_
  end

  # solve ODE for N many particles transiting this factor
  kS = cf.factor.keepSolution
  if !isnothing( kS )
    resize!(kS, N)
    # @info "kS" length(kS)
    for i in 1:N
      kS[i] = (;t=Vector{Float64}(),u=Vector{Vector{Float64}}())
    end
  end
  # solve likely elements
  for i = 1:N
    # TODO, does this respect hyporecipe ???
    idxArr = (k -> cf._legacyParams[k][i]).(1:length(cf._legacyParams))
    oderes = _solveFactorODE!(meas[i], prob, u0pts[i], _maketuplebeyond2args(idxArr...)...)
    if !isnothing( kS )
      resize!(kS[i].t, length(oderes.t))
      resize!(kS[i].u, length(oderes.u))
      for j = 1:length(oderes.t)
        kS[i].t[j] = oderes.t[j]
        kS[i].u[j] = oderes.u[j]
      end
    end
  end

  # return meas, M
  return map(x -> (x, M), meas)
end
# getDimension(oderel.domain)



end # module