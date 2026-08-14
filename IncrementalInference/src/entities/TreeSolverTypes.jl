# TODO: towards solver options and not :parametric symbol
# used for parametric tree solve.

abstract type AbstractTreeSolver end

"""
    $TYPEDEF

Options for the **tangent-space (linearized) parametric** Bayes tree solve.

One sweep is a single Gauss-Newton step about the current reference point, so nonlinear problems need
several; `iters` and `tol` bound that loop.

```julia
solveTree!(fg; algorithm = :parametric)                        # defaults
solveTree!(fg; solver = TangentSpaceSolver(; iters = 25))      # explicit
```

$(TYPEDFIELDS)
"""
Base.@kwdef struct TangentSpaceSolver <: AbstractTreeSolver
  """maximum number of relinearization sweeps"""
  iters::Int = 10
  """convergence threshold on a sweep's largest step"""
  tol::Float64 = 1e-8
  """per-clique re-linearization threshold: a clique reuses its cached Jacobian while none of its
  variables has moved further than this from where the cache was built.  `0` disables the cache and
  re-linearizes every sweep — the baseline the caching must reproduce exactly.

  Where `tol` decides when the *solve* is finished, this decides when a *clique* may skip work.

  This also sets the accuracy floor: Each clique's answer came from a linearization up to 
  `relinearizeTol` away from where it now sits, so the solve converges to roughly `relinearizeTol` 
  rather than `tol` whenever `relinearizeTol > tol`."""
  relinearizeTol::Float64 = 0.0
  """use the exact differential `Jᵣ` when re-anchoring a message onto the receiver's tangent space,
  rather than the flat `A = I`.

  `A = I` is exact on abelian groups and elsewhere wrong at first order in the offset; `Jᵣ` removes that
  term, at one `jacobian_exp` per variable per message.

  Off by default: at `relinearizeTol = 0` every clique shares a reference point, so the offset is zero
  and there is nothing to correct."""
  exactJacobian::Bool = false
end

_algorithmLabel(::TangentSpaceSolver) = :parametric

"""
CSM state transitions one clique runs per parametric sweep.  A rough count used to size the
`limititers` warning and the progress label — not a contract.
"""
const CSM_STEPS_PER_SWEEP = 6

"""
CSM state transitions a clique runs once, outside the sweep loop: recycling check, subgraph build,
presolve checklist, linearized setup, finalize, and the write-back.
"""
const CSM_STEPS_FIXED = 6

"""
    $SIGNATURES

Rough upper bound on total CSM state transitions, for labelling the progress meter only.

The default solve is a single up/down pass, so a fixed per-clique budget covers it.  The parametric
CSM loops *internally* (see `checkConverged_ParametricStateMachine`), so its budget scales with
`solver.iters`.
"""
_approxCSMIters(ncliqs::Int, ::Any) = ncliqs * 24
_approxCSMIters(ncliqs::Int, solver::TangentSpaceSolver) = ncliqs * CSM_STEPS_PER_SWEEP * solver.iters
