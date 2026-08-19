# TODO: towards solver options and not :parametric symbol
# used for parametric tree solve.

abstract type AbstractTreeSolver end

"""
    $TYPEDEF

Options for the tangent-space (linearized) parametric Bayes tree solve.

$(TYPEDFIELDS)
"""
Base.@kwdef struct TangentSpaceSolver <: AbstractTreeSolver
  """maximum number of relinearization sweeps"""
  iters::Int = 10
  """convergence threshold on the largest step across a sweep"""
  tol::Float64 = 1e-8
  """clique reuses its cached Jacobian while no variable has moved further than this from the cache point.
  `0` disables caching. Sets the accuracy floor: converges to `max(tol, relinearizeTol)`."""
  relinearizeTol::Float64 = 0.0
  """use exact `Jᵣ` differential when re-anchoring messages onto the receiver's tangent space.
  Off by default: at `relinearizeTol = 0` anchors coincide and the correction is zero."""
  exactJacobian::Bool = false
end

_algorithmLabel(::TangentSpaceSolver) = :parametric

# rough CSM state transitions per parametric sweep (used to size the `limititers` warning)
const CSM_STEPS_PER_SWEEP = 6

# CSM state transitions run once per clique outside the sweep loop
const CSM_STEPS_FIXED = 6

""" $SIGNATURES
Rough upper bound on total CSM state transitions, for sizing the progress meter.
"""
_approxCSMIters(ncliqs::Int, ::Any) = ncliqs * 24
_approxCSMIters(ncliqs::Int, solver::TangentSpaceSolver) = ncliqs * CSM_STEPS_PER_SWEEP * solver.iters
