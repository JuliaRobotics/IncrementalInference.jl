# ======================================================================================
# Parametric Bayes tree solve — tangent-space (linearized) core.
#   up:   linearize factors -> (J,r); fuse children; Schur-complement frontals; send Λ_S to parent.
#   down: back-substitute ΔF, step x = exp(x0,Δ), pass to children.
# Graph boundary is confined to `_readBasePoint`, `_writeBundle!`, `_cliqueLayout`, and
# `CliqueLinearizer`. Per-sweep functions take no graph.
# Algebra lives in `BeliefAlgebra.jl`.
# ======================================================================================


_getStateTypes(dfg::AbstractDFG, syms::AbstractVector{Symbol}) =
  DataType[typeof(getStateKind(getVariable(dfg, s))) for s in syms]

""" $SIGNATURES
Read variables from `dfg` into a point container ordered by `varlabelsAP`.
"""
function _readBasePoint(subfg::AbstractDFG, varlabelsAP, solveKey::Symbol)
  return map(varlabelsAP) do label
    DFG.refMeans(getState(subfg, label, solveKey))[1]
  end
end

""" $SIGNATURES
Write base points and, for `covarianceLabels` only, covariance blocks from `bundle` back to `dfg`.
Separator covariances are owned by their frontal clique, so passing all labels would be wrong.

TODO fix `refMeans`/`refCovariances`: build `State` objects and use `mergeStates!`.
"""
function _writeBundle!(
  dfg::AbstractDFG,
  bundle::DensityBundlePoint{<:Any, <:TangentNormal},
  solveKey::Symbol;
  covarianceLabels::AbstractVector{Symbol},
)
  Σ = bundle.fibre.Σ
  for v in eachvariable(bundle)
    state = getState(dfg, v.label, solveKey)
    DFG.refMeans(state)[1] = v.point
    v.label in covarianceLabels || continue
    DFG.refCovariances(state)[1] = Σ[v.range, v.range]
  end
  return bundle
end

# ======================================================================================
# Clique elimination
# ======================================================================================

""" $SIGNATURES
Build the clique's `JointLayout`: frontals then separators, each grouped by state type.
Built once per clique; reused across all sweeps.
"""
function _cliqueLayout(
  subfg::AbstractDFG,
  frontals::AbstractVector{Symbol},
  separators::AbstractVector{Symbol},
)
  _, ap, _ = buildPartitionedSolveManifold(
    (getVariable.(subfg, frontals), getVariable.(subfg, separators)),
  )
  labels = Symbol[collect(ap)...]
  return JointLayout(labels, _getStateTypes(subfg, labels), length.(collect(ap.x)))
end

""" $SIGNATURES
Build a [`CliqueLinearizer`](@ref) for `layout`'s variables and `faclabels`.
`p0` sizes internal buffers only; the linearization point is passed at each [`linearize!`](@ref) call.
`layout` must be the clique's own object (`covers` uses pointer equality).
"""
function CliqueLinearizer(
  subfg::AbstractDFG,
  layout::JointLayout,
  faclabels::AbstractVector{Symbol},
  p0 = nothing;
  solveKey::Symbol = :parametric,
  partition = nothing,
)
  M, costF!, jacF! = build_costF_jacF(
    subfg, layout.labels, faclabels; is_sparse = false, solveKey, p0, partition,
  )

  return CliqueLinearizer(
    M, costF!, jacF!,
    zeros(length(jacF!.res), manifold_dimension(M)),
    zeros(length(jacF!.res)),
    layout,
  )
end

""" $SIGNATURES
Linearize the clique's factors, fuse children's upward messages, and Schur-complement the frontals out.
Constructs `cliqdata.elimination` on the first sweep; updates it in place thereafter.
Solve-progress fields (`laststep`, `posterior`) are intentionally not overwritten.
Takes no graph; `p_current` and `childMsgs` are the only per-sweep inputs.
"""
function eliminateCliqueFrontals!(
  cliqdata::BayesTreeNodeData,
  p_current,
  childMsgs::AbstractVector{LinearizedLikelihood} = LinearizedLikelihood[];
  relinearizeTol::Real = 0.0,
  exactJacobian::Bool = false,
)
  layout = cliqdata.cliquelayout
  frontals = cliqdata.frontalIDs
  linearizer = cliqdata.linearizer
  cache = cliqdata.elimination

  @assert covers(linearizer, layout) "linearizer belongs to another clique: $(linearizer.layout.labels) against $(layout.labels)"

  reusedlinearization =
    !isnothing(cache) && _isLinearizationCurrent(cache, p_current, relinearizeTol)

  # pre-fusion system anchored at the Jacobian point; `fuse`/`pullback` are functional so reuse is safe
  cliquelinearization = if reusedlinearization
    cache.cliquelinearization
  elseif isnothing(linearizer)
    # all information arrives through children; start with a zero system
    _n = length(layout)
    DensityBundlePoint(layout, p_current, CotangentNormal(zeros(_n, _n), zeros(_n), 0.0))
  else
    # Gauss-Newton normal equations; `r` is whitened, so logmass = -½‖r‖²
    J, r = linearize!(linearizer, p_current)
    DensityBundlePoint(layout, p_current, CotangentNormal(J' * J, -(J' * r), -0.5 * dot(r, r)))
  end

  # pull each child message back onto this clique's tangent space before fusing
  subtree = cliquelinearization
  for childMsg in childMsgs
    subtree = fuse(subtree, pullback(childMsg.bundle, cliquelinearization, exactJacobian))
  end

  @assert all(subtree.index[s] <= length(frontals) for s in frontals) "clique layout is not frontals-first: $(subtree.labels) against frontals $frontals"

  nfrontal = 0
  for s in frontals
    nfrontal += length(getCoordrange(subtree, s))
  end
  idx_F = 1:nfrontal
  idx_S = (nfrontal + 1):length(subtree)
  seps = subtree.labels[(length(frontals) + 1):end]

  Λ_subtree_mat = subtree.fibre.Λ
  _diagnoseFrontalRank(
    Symmetric(view(Λ_subtree_mat, idx_F, idx_F)),
    view(Λ_subtree_mat, idx_F, idx_S),
    frontals,
  )

  # reuse the separator layout across sweeps; only fibre and base points change
  conditional, separatorlikelihood = eliminate(
    subtree, frontals;
    seps, idx_F, idx_S,
    seplayout = isnothing(cache) ? subsetLayout(layout, seps) : cache.separatorlikelihood.layout,
  )
  @assert conditional.separators == separatorlikelihood.labels

  # `Λ_subtree_mat` aliases `subtree`; `_storeElimination!` copies into its own buffer
  cliqdata.elimination = _storeElimination!(
    cache,
    (;
      idx_F, idx_S, cliquelinearization, Λ_subtree = Λ_subtree_mat,
      separatorlikelihood, conditional, reusedlinearization,
    ),
  )
  return cliqdata.elimination
end

# zero Δ and Λ=0 placeholder; the downward pass overwrites before anything reads it
function _seedPosterior(cliquelinearization::DensityBundlePoint)
  n = length(cliquelinearization)
  return LinearizedBelief(
    DensityBundlePoint(cliquelinearization.layout, cliquelinearization.point, zeros(n, n)),
    zeros(n),
  )
end

# dispatch on Nothing vs CliqueElimination to construct-or-update without branching in the caller
_storeElimination!(::Nothing, e) = CliqueElimination(;
  e.idx_F, e.idx_S, e.cliquelinearization, Λ_subtree = Matrix(e.Λ_subtree),
  e.separatorlikelihood, e.conditional, e.reusedlinearization,
  posterior = _seedPosterior(e.cliquelinearization),
)

function _storeElimination!(elimination::CliqueElimination, e)
  elimination.idx_F = e.idx_F
  elimination.idx_S = e.idx_S
  elimination.cliquelinearization = e.cliquelinearization
  copyto!(elimination.Λ_subtree, e.Λ_subtree)
  elimination.separatorlikelihood = e.separatorlikelihood
  elimination.conditional = e.conditional
  elimination.reusedlinearization = e.reusedlinearization
  return elimination
end

""" $SIGNATURES
Downward pass: back-substitute `ΔS` and recover the exact joint precision via the cavity correction.
Returns `(; ΔF, Δ_full, Λ_exact)` — all anchored at `cliquelinearization.point`; nothing is stepped here.
`parentBelief` is `nothing` at a root (and `ΔS` is empty there).
See `dev/parametric_tree_solve.md §24` for the cavity derivation.
"""
function downsolveClique(
  elimination,
  ΔS::AbstractVector,
  parentBelief::Union{Nothing, DensityBundlePoint} = nothing,
)
  ΔF = backsubstitute(elimination.conditional, ΔS)

  Δ_full = zeros(length(elimination.cliquelinearization))
  Δ_full[elimination.idx_F] .= ΔF
  isempty(ΔS) || (Δ_full[elimination.idx_S] .= ΔS)

  Λ_exact = Matrix(elimination.Λ_subtree)
  if !isnothing(parentBelief)
    Λ_exact[elimination.idx_S, elimination.idx_S] .+= calcCavityprecision(parentBelief, elimination.separatorlikelihood)
  end
  return (; ΔF, Δ_full, Λ_exact = Symmetric(Λ_exact))
end

""" $SIGNATURES
Return `true` if the cached Jacobian at `cache.cliquelinearization.point` is still valid at `p_current`.
`relinearizeTol = 0` forces re-linearization every sweep.
"""
function _isLinearizationCurrent(cache::CliqueElimination, p_current, relinearizeTol::Real)
  relinearizeTol <= 0 && return false
  bundle = cache.cliquelinearization
  layout = bundle.layout
  length(p_current.x) == length(layout.blocks) || return false
  for b in eachindex(layout.blocks)
    G = getManifold(layout.blocks[b].statetype)
    ref, cur = bundle.point.x[b], p_current.x[b]
    length(cur) == length(ref) || return false
    for i in eachindex(cur)
      norm(log_coord(G, ref[i], cur[i])) > relinearizeTol && return false
    end
  end
  return true
end

""" $SIGNATURES
Re-express tangent deltas of `syms` from `msg` at `own`'s reference points:
`ΔS_local = log(p_own, exp(p_msg, Δ_msg))`.
Short-circuits when both anchors coincide (measured ~68% of calls at `relinearizeTol=0`).
"""
function _reanchorDeltas(
  msg::LinearizedContent,
  syms::AbstractVector{Symbol},
  own::DensityBundlePoint,
)
  msgBelief = msg.bundle  # covers variables the local subgraph may lack
  for s in syms
    @assert hasLabel(msgBelief, s) "downward message has no delta for $s"
  end

  n = 0
  for s in syms
    n += length(getCoordrange(msgBelief, s))
  end
  out = Vector{Float64}(undef, n)
  basis = DefaultLieAlgebraOrthogonalBasis()
  offset = 0
  for s in syms
    r = getCoordrange(msgBelief, s)
    Δ_msg = view(msg.Δ, r)
    target = view(out, (offset + 1):(offset + length(r)))
    offset += length(r)

    if !hasLabel(own, s) || getBasepoint(msgBelief, s) == getBasepoint(own, s)
      target .= Δ_msg
      continue
    end

    G = getManifold(getStatetype(msgBelief.layout, msgBelief.index[s]))
    target .= log_coord(
      G,
      getBasepoint(own, s),
      exp_coord(G, getBasepoint(msgBelief, s), Δ_msg, basis),
      basis
    )
  end
  return out
end

# ======================================================================================
# Outer relinearization loop
# ======================================================================================

""" $SIGNATURES
Run the tangent-space parametric tree solve; returns `(tree, sweeps, status)`.
Sweeps Gauss-Newton steps until convergence or `solver.iters` is reached.

!!! warning "Reuse is only valid while the elimination structure holds"
    A reused `tree` must still match `fg`; adding or removing variables or factors invalidates it.
    Pass `nothing` (the default) to rebuild.
"""
function solveTreeParametric!(
  fg::AbstractDFG,
  tree::Union{Nothing, AbstractBayesTree} = nothing;
  solver::TangentSpaceSolver = TangentSpaceSolver(),
  solveKey::Symbol = :parametric,
  eliminationOrder::Union{Nothing, Vector{Symbol}} = nothing,
  eliminationConstraints::Vector{Symbol} = Symbol[],
  kwargs...,
)
  ensureSolvable!(fg)

  # FIXME `limititers` caps state transitions, not sweeps; too low deadlocks the tree
  _needed = CSM_STEPS_PER_SWEEP * solver.iters + CSM_STEPS_FIXED
  if 0 < getSolverParams(fg).limititers < _needed
    @error "FIXME `limititers = $(getSolverParams(fg).limititers)` is below the ~$_needed state transitions \
`solver.iters = $(solver.iters)` sweeps need; the solve might hang."
  end

  makeSolverData!(fg; solveKey, parametric = true)
  if getSolverParams(fg).graphinit
    @info "Ensure variables are all initialized (parametric graphinit on $solveKey)"
    autoinitParametric!(fg; solveKey)
  end

  tree_ = @something(tree, buildSolveTree!(fg; eliminationOrder, eliminationConstraints))

  # reset laststep so a reused tree doesn't appear converged from a prior call
  for (_, cliq) in getCliques(tree_)
    cliqueCache = getCliqueData(cliq).elimination
    isnothing(cliqueCache) || (cliqueCache.laststep = Inf)
  end

  solveTreePass!(fg, tree_; solver, algorithm = :parametric, solveKey, kwargs...)

  sweeps, status = _treeSolveOutcome(tree_)
  @info "solveTreeParametric! finished after $sweeps sweeps, $status"
  return tree_, sweeps, status
end

""" $SIGNATURES
Return `(sweeps, status)` of the last parametric solve from the tree's roots.
"""
function _treeSolveOutcome(tree::AbstractBayesTree)
  sweeps = 0
  status = CONVERGED
  for (_, cliq) in getCliques(tree)
    isempty(getParent(tree, cliq)) || continue
    sweeps = max(sweeps, getCliqueData(cliq).parIter)
    # ERROR_STATUS outranks ITERLIMIT
    cliqstatus = getCliqueStatus(cliq)
    cliqstatus === ITERLIMIT && status !== ERROR_STATUS && (status = ITERLIMIT)
    cliqstatus === ERROR_STATUS && (status = ERROR_STATUS)
  end
  return sweeps, status
end
