# ======================================================================================
# Clique state machine functions — parametric (linearized) path.
# ======================================================================================

"""
    $SIGNATURES

Build the clique's residual/Jacobian machinery — CSM function, run **once** per state machine.
"""
function setupLinear_ParametricStateMachine(csmc::CliqStateMachineContainer)
  logCSM(csmc, "Par-0c $(csmc.cliq.id): setting up linearized solve")

  subfg = csmc.cliqSubFg
  frontals = getCliqFrontalVarIds(csmc.cliq)
  separators = getCliqSeparatorVarIds(csmc.cliq)
  factorLabels = getCliqFactorIdsAll(csmc.cliq)
  cliqdata = getCliqueData(csmc.cliq)

  cliqdata.cliquelayout = @something(cliqdata.cliquelayout, _cliqueLayout(subfg, frontals, separators))

  cliqdata.linearizer = if isempty(factorLabels)
    nothing
  else
    CliqueLinearizer(
      subfg, cliqdata.cliquelayout, factorLabels;
      solveKey = csmc.solveKey, partition = (frontals, separators),
    )
  end

  # re-seed from the subgraph, so an edit made from outside between solves is not hidden by a
  # carried point
  if !isnothing(cliqdata.elimination)
    own = cliqdata.elimination.cliquelinearization
    n = length(own)
    cliqdata.elimination.posterior = LinearizedBelief(
      DensityBundlePoint(
        own.layout,
        _readBasePoint(subfg, getLabelpartition(own), csmc.solveKey),
        zeros(n, n),
      ),
      zeros(n),
    )
  end

  return waitForUp_StateMachine
end

"""
    $SIGNATURES

Upward pass — CSM function.  Linearize, fuse the children's messages, eliminate the frontals, send the
separator system up.  A root has no separators and solves on the way down.
"""
function solveUp_ParametricStateMachine(csmc::CliqStateMachineContainer)
  infocsm(csmc, "Par-3, Solving Up")
  setCliqueDrawColor!(csmc.cliq, "red")

  subfg = csmc.cliqSubFg
  separators = getCliqSeparatorVarIds(csmc.cliq)

  childMsgs = LinearizedLikelihood[]
  for (_, upmsg) in getMessageBuffer(csmc.cliq).upRx
    isnothing(upmsg.linearized) && continue
    push!(childMsgs, upmsg.linearized)
  end

  relinearizeTol = csmc.solver isa TangentSpaceSolver ? csmc.solver.relinearizeTol : 0.0
  exactJacobian = csmc.solver isa TangentSpaceSolver ? csmc.solver.exactJacobian : false
  tol = csmc.solver isa TangentSpaceSolver ? csmc.solver.tol : 0.0
  cliqdata = getCliqueData(csmc.cliq)
  cache = cliqdata.elimination
  childRelinearized = any(childMsg -> childMsg.relinearized, childMsgs)
  childMoving = any(childMsg -> childMsg.moving, childMsgs)

  # the clique's own estimate: carried between sweeps, read back from storage only on the first
  p_current = isnothing(cache) ?
    _readBasePoint(subfg, getLabelpartition(cliqdata.cliquelayout), csmc.solveKey) :
    stepBasepoint(cache.posterior.bundle, cache.posterior.Δ)
  linearizationCurrent =
    !isnothing(cache) && _isLinearizationCurrent(cache, p_current, relinearizeTol)

  elimination, relinearized = if linearizationCurrent && !childRelinearized
    logCSM(csmc, "$(csmc.cliq.id): up skipped, cached factorization still valid")
    (cache, false)
  else
    # constructs on the sweep that has none, updates in place after — and the update leaves `laststep`
    # alone, which is what keeps PROGRESS from being reset by a re-factorization
    newElimination =
      eliminateCliqueFrontals!(cliqdata, p_current, childMsgs; relinearizeTol, exactJacobian)
    logCSM(
      csmc,
      "$(csmc.cliq.id): up eliminated, $(length(childMsgs)) child systems, jacobian $(newElimination.reusedlinearization ? "reused" : "recomputed")",
    )
    (newElimination, true)
  end

  # must describe THIS pass, not whichever pass last wrote the cache
  elimination.relinearized = relinearized
  moving = elimination.laststep > tol || childMoving
  elimination.moving = moving

  _dbgCSMSaveSubFG(csmc, "fg_beforeupsolve")

  if !isempty(separators)
    beliefMsg = LikelihoodMessage(;
      sender = (; id = csmc.cliq.id.value, step = csmc._csm_iter),
      status = UPSOLVED,
      variableOrder = separators,
      msgType = LinearizedMessage(),
      linearized = LinearizedLikelihood(elimination.separatorlikelihood, relinearized, moving),
    )

    getMessageBuffer(csmc.cliq).upTx = beliefMsg
    for edge in getEdgesParent(csmc.tree, csmc.cliq)
      logCSM(csmc, "$(csmc.cliq.id): put! on edge $(edge)")
      putBeliefMessageUp!(csmc.tree, edge, beliefMsg)
    end
  end

  return waitForDown_StateMachine
end

"""
    $SIGNATURES
The parent's belief and this clique's separator deltas, re-expressed onto this clique's reference
points.  A root has neither.
"""
function _downwardInputs(csmc::CliqStateMachineContainer, elimination::CliqueElimination)
  length(getParent(csmc.tree, csmc.cliq)) == 0 && return (nothing, Float64[])
  downmsg = getMessageBuffer(csmc.cliq).downRx
  @assert !isnothing(downmsg) && !isnothing(downmsg.linearized) "no linearized downward message in clique $(csmc.cliq.id.value)"
  # NOTE separator order is the PARTITION's, which `idx_S`, `W` and `separatorlikelihood.fibre.Λ`
  # were all built in
  sepLabels = elimination.separatorlikelihood.labels
  return (
    downmsg.linearized.bundle,
    _reanchorDeltas(downmsg.linearized, sepLabels, elimination.cliquelinearization),
  )
end

"""
    $SIGNATURES

Downward pass — CSM function.  Back-substitute against the parent's deltas, forward the result to the
children, and broadcast the root's termination verdict.  The mathematics is in
[`downsolveClique`](@ref).
"""
function solveDown_ParametricStateMachine(csmc::CliqStateMachineContainer)
  infocsm(csmc, "Lin-5, Solving down (linearized)")
  setCliqueDrawColor!(csmc.cliq, "red")

  subfg = csmc.cliqSubFg
  isroot = length(getParent(csmc.tree, csmc.cliq)) == 0

  # reuse the upward pass's factorization — recomputing would double the clique cost, Jacobian included
  cliqdata = getCliqueData(csmc.cliq)
  cache = cliqdata.elimination
  elimination = if isnothing(cache)
    @warn "linearized down solve: no cached factorization for clique $(csmc.cliq.id.value), recomputing"
    childMsgs = LinearizedLikelihood[]
    for (_, upmsg) in getMessageBuffer(csmc.cliq).upRx
      isnothing(upmsg.linearized) && continue
      push!(childMsgs, upmsg.linearized)
    end
    eliminateCliqueFrontals!(
      cliqdata,
      _readBasePoint(subfg, getLabelpartition(cliqdata.cliquelayout), csmc.solveKey), childMsgs;
      exactJacobian = csmc.solver isa TangentSpaceSolver ? csmc.solver.exactJacobian : false,
    )
  else
    cache
  end

  parentBelief, ΔS = _downwardInputs(csmc, elimination)

  # Decided at the ROOT and broadcast, never locally — a clique that stopped on its own judgement
  # would leave a sibling blocked forever.
  status = if isroot
    csmc.parIter += 1
    getCliqueData(csmc.cliq).parIter = csmc.parIter
    iters = csmc.solver isa TangentSpaceSolver ? csmc.solver.iters : 1
    tol = csmc.solver isa TangentSpaceSolver ? csmc.solver.tol : 0.0
    relin = getCliqueData(csmc.cliq).elimination.relinearized
    still_moving = getCliqueData(csmc.cliq).elimination.moving
    if !still_moving || !relin
      why = !relin ? "nothing re-linearized" : "step < tol"
      logCSM(csmc, "$(csmc.cliq.id): CONVERGED after $(csmc.parIter) sweeps ($why)")
      CONVERGED
    elseif csmc.parIter >= iters
      # NOT `CONVERGED`: it stopped because it ran out, and a log should not make you infer which
      logCSM(csmc, "$(csmc.cliq.id): ITERLIMIT at $(csmc.parIter) sweeps, still re-linearizing")
      ITERLIMIT
    else
      DOWNSOLVED
    end
  else
    getMessageBuffer(csmc.cliq).downRx.status
  end
  setCliqueStatus!(csmc.cliq, status)

  down = downsolveClique(elimination, ΔS, parentBelief)
  getCliqueData(csmc.cliq).elimination.laststep = norm(down.ΔF)

  _dbgCSMSaveSubFG(csmc, "fg_afterdownsolve_$(csmc.parIter)")

  # kept rather than dropped: the children receive it and `finalizeLinear` converts it to the
  # covariance without re-deriving it.  `η === nothing` — the mean travels as `Δ_full`, which survives
  # a singular `Λ_exact` where `η = ΛΔ` would silently drop the minimum-norm choice.
  posterior = LinearizedBelief(
    DensityBundlePoint(
      elimination.cliquelinearization.layout, elimination.cliquelinearization.point,
      Matrix(down.Λ_exact),
    ),
    down.Δ_full,
  )
  getCliqueData(csmc.cliq).elimination.posterior = posterior

  beliefMsg = LikelihoodMessage(;
    sender = (; id = csmc.cliq.id.value, step = csmc._csm_iter),
    status = status,
    variableOrder = elimination.cliquelinearization.labels,
    msgType = LinearizedMessage(),
    linearized = posterior,
  )

  @sync for edge in getEdgesChildren(csmc.tree, csmc.cliq)
    logCSM(csmc, "$(csmc.cliq.id): put! on edge $(edge)")
    @async putBeliefMessageDown!(csmc.tree, edge, beliefMsg)
  end

  logCSM(csmc, "$(csmc.cliq.id): linearized sweep completed, status $status")

  return checkConverged_ParametricStateMachine
end

"""
    $SIGNATURES
Loop or finish — CSM function closing one sweep.  Terminates only on a status the root broadcast.
"""
function checkConverged_ParametricStateMachine(csmc::CliqStateMachineContainer)
  status = getCliqueStatus(csmc.cliq)
  if status in (CONVERGED, ITERLIMIT)
    setCliqueDrawColor!(csmc.cliq, status === CONVERGED ? "green" : "orange")
    logCSM(csmc, "$(csmc.cliq.id): finishing on $status")
    return finalizeLinear_ParametricStateMachine
  end
  logCSM(csmc, "$(csmc.cliq.id): another sweep")
  return waitForUp_StateMachine
end

"""
    $SIGNATURES

Write the converged clique back into its subgraph — CSM function, run **once** per state machine.
Solves nothing: the last downward pass left the answer in `elimination.posterior`.
"""
function finalizeLinear_ParametricStateMachine(csmc::CliqStateMachineContainer)
  elimination = getCliqueData(csmc.cliq).elimination
  posterior = isnothing(elimination) ? nothing : elimination.posterior
  if isnothing(posterior)
    logCSM(csmc, "$(csmc.cliq.id): nothing to finalize, no downward solve ran")
    return updateFromSubgraph_StateMachine
  end

  # convert BEFORE stepping: a covariance must never be made by inverting a transported precision.
  # Frontals only, but it takes the whole joint to produce them — `Σ_FF` is a block of `Λ⁻¹`.
  _writeBundle!(
    csmc.cliqSubFg,
    stepBelief(calcMomentform(posterior.bundle), posterior.Δ),
    csmc.solveKey;
    covarianceLabels = getCliqFrontalVarIds(csmc.cliq),
  )

  return updateFromSubgraph_StateMachine
end
