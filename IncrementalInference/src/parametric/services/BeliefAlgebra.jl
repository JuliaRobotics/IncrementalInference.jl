# ======================================================================================
# Belief algebra — the operations that make up a tree solve, on bundles and fibres alone.
# Nothing here takes a graph, a message, or a clique.
# ======================================================================================

"""
    $SIGNATURES

Inverse of a symmetric PSD matrix on its observable subspace.
"""
_invPSD(A::AbstractMatrix; rtol::Real = 1e-9) = _factorPSD(A; rtol).inv

""" $SIGNATURES
Apply a PSD factor's inverse to `B`: triangular solve at full rank, multiply at singular.
"""
_applyInverse(C::CholeskyPivoted, B::AbstractVecOrMat) = C \ B
_applyInverse(Minv::AbstractMatrix, B::AbstractVecOrMat) = Minv * B

""" $SIGNATURES
Pivoted Cholesky rank test; returns `(rank, factor)`.

!!! warning "`tol` is not optional"
    Scale `rtol` by the largest diagonal to make this a rank test, not just a positive-pivot check.
"""
function _rankPSD(M::Symmetric; rtol::Real = 1e-9)
  n = size(M, 1)
  n == 0 && return (rank = 0, factor = nothing)
  scale = maximum(diag(M))
  scale > 0 || return (rank = 0, factor = nothing)
  C = cholesky(M, RowMaximum(); tol = rtol * scale, check = false)
  return (rank = C.rank, factor = C)
end

""" $SIGNATURES
Log pseudo-determinant and rank of a symmetric PSD matrix on its observable subspace.
"""
function _logDetPSD(A::AbstractMatrix; rtol::Real = 1e-9)
  M = Symmetric(A)
  n = size(M, 1)
  n == 0 && return (0.0, 0)

  r, C = _rankPSD(M; rtol)
  r == n && return (2 * sum(log, diag(C.U)), n)

  λ = eigvals(M)
  λmax = maximum(abs, λ)
  keep = filter(x -> x > rtol * max(λmax, eps()), λ)
  return (isempty(keep) ? 0.0 : sum(log, keep), length(keep))
end

""" $SIGNATURES
Inverse, log pseudo-determinant, and rank of a symmetric PSD matrix from one decomposition.
Avoids decomposing the same block twice when both inverse and determinant are needed.
"""
function _factorPSD(A::AbstractMatrix; rtol::Real = 1e-9)
  M = Symmetric(A)
  n = size(M, 1)
  n == 0 && return (inv = Symmetric(zeros(0, 0)), logdet = 0.0, rank = 0)

  r, C = _rankPSD(M; rtol)
  r == n && return (inv = Symmetric(inv(C)), logdet = 2 * sum(log, diag(C.U)), rank = n)

  # singular: one eigen serves both pseudo-inverse and pseudo-determinant on the same subspace
  eig = eigen(M)
  λmax = maximum(abs, eig.values)
  keep = eig.values .> rtol * max(λmax, eps())
  λ = eig.values[keep]
  isempty(λ) && return (inv = Symmetric(zeros(n, n)), logdet = 0.0, rank = 0)
  V = eig.vectors[:, keep]
  return (inv = Symmetric(V * Diagonal(1 ./ λ) * V'), logdet = sum(log, λ), rank = length(λ))
end

## ======================================================================================
## Geometry over a bundle's layout
## ======================================================================================

""" $SIGNATURES
`exp` at `p` from tangent coordinates `c`. Uses `get_vector` (not `hat`) for point-aware tangent
construction; explicit basis avoids silent coordinate disagreements with `hat`/`vee`.
"""
exp_coord(G, p, c, B::AbstractBasis = DefaultLieAlgebraOrthogonalBasis()) =
  exp(G, p, get_vector(G, p, c, B))

""" $SIGNATURES
In-place [`exp_coord`](@ref) writing into mutable `q`.
"""
exp_coord!(G, q, p, c, B::AbstractBasis = DefaultLieAlgebraOrthogonalBasis()) =
  exp!(G, q, p, get_vector(G, p, c, B))

""" $SIGNATURES
Coordinates of `log(p, q)`. Inverse of [`exp_coord`](@ref).
"""
log_coord(G, p, q, B::AbstractBasis = DefaultLieAlgebraOrthogonalBasis()) =
  get_coordinates(G, p, log(G, p, q), B)

""" $SIGNATURES
Block-diagonal `Jᵣ(d)` over `bundle`'s layout — the pullback matrix `A`. Not `Jᵣ(d)⁻¹`.
"""
function LieGroups.jacobian_exp(
  bundle::DensityBundlePoint,
  d::AbstractVector,
  B = DefaultLieAlgebraOrthogonalBasis(),
)
  length(d) == length(bundle) ||
    error("jacobian_exp: d has length $(length(d)) but $(bundle.labels) span $(length(bundle)) dimensions")
  layout = bundle.layout
  J = Matrix(1.0I, length(bundle), length(bundle))
  for b in eachindex(layout.blocks)
    G = getManifold(layout.blocks[b].statetype)
    start = layout.blockstarts[b]
    block = bundle.point.x[b]
    for i in eachindex(block)
      r = layout.ranges[start + i - 1]
      all(iszero, view(d, r)) && continue   # Jᵣ(0) = I
      J[r, r] .= jacobian_exp(G, get_vector(G, block[i], view(d, r), B), B)
    end
  end
  return J
end

""" $SIGNATURES
Bundle point over `syms` with `bundle`'s base points and the given fibre.
`syms` must be an order-preserving subset of `bundle.labels`. Pass `layout` to avoid recomputing it.
"""
function subsetBundle(
  bundle::DensityBundlePoint,
  syms::AbstractVector{Symbol},
  fibre;
  layout::JointLayout = subsetLayout(bundle.layout, syms),
)
  keep = Set(syms)
  own = bundle.layout
  blocks = map(eachindex(own.blocks)) do b
    start = own.blockstarts[b]
    block = bundle.point.x[b]
    [block[i] for i in eachindex(block) if own.labels[start + i - 1] in keep]
  end
  return DensityBundlePoint(layout, ArrayPartition(filter(!isempty, blocks)...), fibre)
end

""" $SIGNATURES
Positive-definite covariance from `fibre`, or `nothing` if the precision is not positive-definite.
"""
function _properCovariance(fibre::CotangentNormal)
  length(fibre) == 0 && return zeros(0, 0)
  F = cholesky(Symmetric(fibre.Λ); check = false)
  issuccess(F) || return nothing
  return inv(F)
end

function _properCovariance(fibre::TangentNormal)
  length(fibre) == 0 && return zeros(0, 0)
  issuccess(cholesky(Symmetric(fibre.Σ); check = false)) || return nothing
  return fibre.Σ
end

""" $SIGNATURES
`true` if the tangent-space Gaussian is a sound approximation on the manifold (injectivity + volume
tests at `nsigma` σ). Not yet called.
"""
function isConcentrated(
  bundle::DensityBundlePoint;
  nsigma::Real = 3,
  tol::Real = 1e-2,
  atol::Real = 1e-6,
)
  Σ = _properCovariance(bundle.fibre)
  isnothing(Σ) && return false
  B = DefaultLieAlgebraOrthogonalBasis()
  for v in eachvariable(bundle)
    isempty(v.range) && continue
    eig = eigen(Symmetric(Σ[v.range, v.range]))
    G = getManifold(v.statetype)
    𝔤 = LieAlgebra(G)
    p = v.point
    for i in eachindex(eig.values)
      eig.values[i] > 0 || return false
      coords = (nsigma * sqrt(eig.values[i])) .* eig.vectors[:, i]
      isapprox(log_coord(G, p, exp_coord(G, p, coords, B), B), coords; atol) || return false
      abs(det(jacobian_exp(G, get_vector(G, p, coords, B), B)) - 1) <= tol || return false
    end
  end
  return true
end

## ======================================================================================
## Moving a belief between base points
## ======================================================================================

""" $SIGNATURES
Stacked log-coordinates from `own`'s base points to `bundle`'s: `d[v] = log(p_own, p_bundle)`.
Zero for labels `own` does not hold (treated as co-anchored).
"""
function logCoordinates(bundle::DensityBundlePoint, own::DensityBundlePoint)
  d = zeros(length(bundle))
  for v in eachvariable(bundle)
    hasLabel(own, v.label) || continue
    G = getManifold(v.statetype)
    d[v.range] .= log_coord(G, getBasepoint(own, v.label), v.point)
  end
  return d
end

""" $SIGNATURES
Pull back a [`CotangentNormal`](@ref) by offset `d` and Jacobian `A` (default `I` = flat pullback):
`Λ' = AᵀΛA`, `η' = Aᵀ(η + Λd)`, `logmass' = logmass − ηᵀd − ½dᵀΛd`.
"""
function pullback(fibre::CotangentNormal, d::AbstractVector, A = I)
  iszero(d) && A === I && return fibre
  Λ = fibre.Λ 
  η = getInfovector(fibre)
  logmass_prime = fibre.logmass - dot(η, d) - 0.5 * dot(d, Λ * d)
  A === I && return CotangentNormal(Λ, η + Λ * d, logmass_prime)
  return CotangentNormal(A' * Λ * A, A' * (η + Λ * d), logmass_prime)
end

""" $SIGNATURES
Pull `bundle` back onto `own`'s tangent space, deriving offset and Jacobian internally.
`exactJacobian = true` uses `A = Jᵣ(d)`; default is the flat pullback `A = I`.
"""
function pullback(
  bundle::DensityBundlePoint,
  own::DensityBundlePoint{PT},
  exactJacobian::Bool = false,
) where PT
  d = logCoordinates(bundle, own)
  A = exactJacobian ? jacobian_exp(bundle, d) : I
  layout = bundle.layout
  blocks = map(Tuple(eachindex(layout.blocks))) do b
    start = layout.blockstarts[b]
    block = bundle.point.x[b]
    map(eachindex(block)) do i
      s = layout.labels[start + i - 1]
      hasLabel(own, s) ? getBasepoint(own, s) : block[i]
    end
  end
  return DensityBundlePoint(bundle.layout, ArrayPartition(blocks), pullback(bundle.fibre, d, A))
end

""" $SIGNATURES
Step `bundle`'s base points by `Δ` and transport the fibre under `Jᵣ(Δ)⁻¹`.
`bundle` must be centred.
"""
function stepBelief(bundle::DensityBundlePoint, Δ::AbstractVector)
  isCentred(bundle) ||
    error("stepBelief: fibre must be centred — the mean is the `Δ` argument, not `η`")
  return DensityBundlePoint(
    bundle.layout, stepBasepoint(bundle, Δ),
    _transportFibre(bundle.fibre, jacobian_exp(bundle, Δ)),
  )
end

""" $SIGNATURES
Step `bundle`'s base points by `Δ`, returning only the new point container.
Use [`stepBelief`](@ref) when the fibre must be carried.
"""
function stepBasepoint(bundle::DensityBundlePoint, Δ::AbstractVector)
  length(Δ) == length(bundle) ||
    error("stepBasepoint: Δ has length $(length(Δ)) but $(bundle.labels) span $(length(bundle)) dimensions")
  layout = bundle.layout
  blocks = map(eachindex(layout.blocks)) do b
    _stepBlock(layout, bundle.point.x[b], Δ, b)
  end
  return ArrayPartition(blocks...)
end

""" $SIGNATURES
Step block `b`'s base points by their slice of `Δ`. `G` resolves once per block (not per variable).
Uses `get_vector`/`exp!` with one mutable buffer `q` reused across variables to avoid allocation.
"""
function _stepBlock(layout::JointLayout, block::AbstractVector, Δ::AbstractVector, b::Int)
  isempty(block) && return similar(block)
  G = getManifold(layout.blocks[b].statetype)
  basis = DefaultLieAlgebraOrthogonalBasis()
  start = layout.blockstarts[b]
  out = similar(block)
  q = allocate(block[1])  # mutable buffer; `exp!` writes here, then we convert to immutable
  for i in eachindex(block)
    p = block[i]
    exp_coord!(G, q, p, view(Δ, layout.ranges[start + i - 1]), basis)
    out[i] = convert(eltype(block), q)
  end
  return out
end

""" $SIGNATURES
Transport a fibre across a base-point step with Jacobian `J = Jᵣ(Δ)`.
Covariance: `JΣJᵀ`; precision: `AᵀΛA` with `A = J⁻¹`. Dispatch prevents mixing them.
"""
_transportFibre(fibre::CotangentNormal, J::AbstractMatrix) =
  (A = inv(J); CotangentNormal(A' * fibre.Λ * A, nothing, fibre.logmass))

_transportFibre(fibre::TangentNormal, J::AbstractMatrix) =
  TangentNormal(J * fibre.Σ * J', nothing, fibre.logmass)

## ======================================================================================
## Combining and reducing
## ======================================================================================

""" $SIGNATURES
Scatter-add `source`'s fibre into `target`'s by label. `source` must be pulled back onto `target`'s
tangent space first; alignment is by label, so ordering does not need to match.
"""
function fuse(target::DensityBundlePoint, source::DensityBundlePoint)
  idx = getCoordindices(target, source.labels)
  length(idx) == length(source) ||
    error("fuse: source spans $(length(source)) coordinates but maps to $(length(idx)) in target")
  Λ = copy(target.fibre.Λ)
  η = copy(getInfovector(target))
  Λ[idx, idx] .+= source.fibre.Λ
  η[idx] .+= getInfovector(source)
  return DensityBundlePoint(
    target.layout,
    target.point,
    CotangentNormal(Λ, η, target.fibre.logmass + source.fibre.logmass),
  )
end

""" $SIGNATURES
Marginalize `bundle`'s fibre to `syms` via Schur complement. `syms` is what survives.
"""
function marginalizeTo(bundle::DensityBundlePoint, syms::AbstractVector{Symbol})
  return marginalizeTo(bundle.fibre, getCoordindices(bundle, syms))
end

""" $SIGNATURES
Schur-complement marginalization onto coordinate indices `keep`. Returns the same fibre form.
"""
function marginalizeTo(fibre::CotangentNormal, keep::AbstractVector{Int})
  Λ = fibre.Λ
  drop = setdiff(1:length(fibre), keep)
  isempty(drop) && return CotangentNormal(
    Λ[keep, keep], isnothing(fibre.η) ? nothing : fibre.η[keep], fibre.logmass,
  )

  Λ_dd = Λ[drop, drop]
  Λ_dd_inv, logdet_dd, rank_dd = _factorPSD(Λ_dd)
  Λ_keep = Symmetric(Λ[keep, keep] - Λ[keep, drop] * Λ_dd_inv * Λ[drop, keep])

  η_keep = isnothing(fibre.η) ? nothing : fibre.η[keep] - Λ[keep, drop] * (Λ_dd_inv * fibre.η[drop])

  η_drop = isnothing(fibre.η) ? nothing : fibre.η[drop]
  logmass = fibre.logmass - 0.5 * logdet_dd + 0.5 * rank_dd * log(2π)
  isnothing(η_drop) || (logmass += 0.5 * dot(η_drop, Λ_dd_inv * η_drop))

  return CotangentNormal(Λ_keep, η_keep, logmass)
end

""" $SIGNATURES
Convert canonical fibre to moment form. Unobservable directions get zero variance (truth is infinite).
"""
function calcMomentform(bundle::DensityBundlePoint{<:Any, <:CotangentNormal})
  Σ = _invPSD(bundle.fibre.Λ)
  return DensityBundlePoint(
    bundle.layout, bundle.point,
    TangentNormal(Σ, isCentred(bundle) ? nothing : Σ * bundle.fibre.η, bundle.fibre.logmass),
  )
end

## ======================================================================================
## Factorization: p(F,S) = p(F│S)·p(S)
## ======================================================================================

""" $SIGNATURES
Diagnose a rank-deficient frontal block: `:full`, `:benign`, or `:malignant`.
Benign = null directions decoupled from separators; malignant = upward message would over-claim.
"""
function _diagnoseFrontalRank(
  Λ_FF::AbstractMatrix,
  Λ_FS::AbstractMatrix,
  ctx = "";
  rtol::Real = 1e-9,
  coupling_rtol::Real = 1e-6,
)
  n = size(Λ_FF, 1)
  M_FF = Symmetric(Λ_FF)

  _rankPSD(M_FF; rtol).rank == n && return :full

  eig = eigen(M_FF)
  λmax = isempty(eig.values) ? 0.0 : maximum(abs, eig.values)
  nullidx = findall(v -> v <= rtol * max(λmax, eps()), eig.values)
  isempty(nullidx) && return :full

  nnull = length(nullidx)
  scale = isempty(Λ_FS) ? 0.0 : opnorm(Λ_FS)
  rel = if isempty(Λ_FS) || scale <= eps()
    0.0
  else
    maximum(k -> norm(Λ_FS' * view(eig.vectors, :, nullidx[k])), eachindex(nullidx)) / scale
  end

  mostlyUndetermined = 2 * nnull > n
  if mostlyUndetermined || rel > coupling_rtol
    why = mostlyUndetermined ? "rank $(n - nnull) of $n" :
          "coupling $(round(rel; sigdigits = 3))"
    @error """clique $(ctx): frontals undetermined ($why) — the eliminated upward message will \
over-claim information.  This clique cannot be summarized in isolation: merge it into its parent, or \
use `eliminationConstraints` so the relevant variables share a clique.""" maxlog = 5
    return :malignant
  end
  @debug "clique $(ctx): $nnull of $n frontal directions undetermined but decoupled from the separators — pseudo-inverse elimination is exact"
  return :benign
end

""" $SIGNATURES
Eliminate `frontals` from `bundle`, returning `(conditional, marginal)` — the factorization `p(F,S) = p(F|S)·p(S)`.
Marginal goes up to the parent; conditional completes on the way down ([`backsubstitute`](@ref)).
Pass `seps`, `idx_F`, `idx_S`, `seplayout` to avoid recomputing structure across sweeps.
"""
function eliminate(
  bundle::DensityBundlePoint,
  frontals::AbstractVector{Symbol};
  seps::AbstractVector{Symbol} = filter(s -> !(s in frontals), bundle.labels),
  idx_F = getCoordspan(bundle, frontals),
  idx_S = getCoordspan(bundle, seps),
  seplayout::JointLayout = subsetLayout(bundle.layout, seps),
)

  Λ, η = bundle.fibre.Λ, getInfovector(bundle)
  Λ_FF = Symmetric(view(Λ, idx_F, idx_F))
  Λ_FS = view(Λ, idx_F, idx_S)
  η_F = view(η, idx_F)

  n_F = size(Λ_FF, 1)
  r, C = _rankPSD(Λ_FF)
  Λ_FF_inv, logdet_FF, rank_FF = if r == n_F
    (C, 2 * sum(log, diag(C.U)), n_F)
  else
    F = _factorPSD(Λ_FF)
    (F.inv, F.logdet, F.rank)
  end
  W = _applyInverse(Λ_FF_inv, Λ_FS)
  v = _applyInverse(Λ_FF_inv, η_F)

  logmass_S = bundle.fibre.logmass + 0.5 * dot(η_F, v) - 0.5 * logdet_FF + 0.5 * rank_FF * log(2π)

  marginal = subsetBundle(
    bundle, seps,
    CotangentNormal(
      Symmetric(view(Λ, idx_S, idx_S) - Λ_FS' * W), view(η, idx_S) - Λ_FS' * v, logmass_S,
    );
    layout = seplayout,
  )
  conditional = GaussianConditional(Λ_FF, v, W, collect(Symbol, frontals), seps)
  return (conditional, marginal)
end

"""
    $SIGNATURES
Complete the conditional against solved separator deltas: `ΔF = v - W·ΔS`.

`ΔS` must be ordered by `c.separators`.  Empty `ΔS` is a root, where the conditional is already the
answer.
"""
backsubstitute(c::GaussianConditional, ΔS::AbstractVector) =
  isempty(ΔS) ? copy(c.v) : c.v - c.W * ΔS

""" $SIGNATURES
Precision outside this subtree over the separators: `Λ_cavity = Λ_S^parent − Λ_S^up`.
A subtraction: the parent marginal already contains this clique's upward message.
"""
calcCavityprecision(parentBelief::DensityBundlePoint, separatorlikelihood::DensityBundlePoint) =
  marginalizeTo(parentBelief, separatorlikelihood.labels).Λ - separatorlikelihood.fibre.Λ

## ======================================================================================
## Evaluation
## ======================================================================================

""" $SIGNATURES
Negative log density `φ(ξ) = ½ξᵀΛξ − ηᵀξ − logmass` at tangent coordinates `ξ`.
"""
function energy(fibre::CotangentNormal, ξ::AbstractVector)
  length(ξ) == length(fibre) ||
    error("energy: ξ has length $(length(ξ)) but the fibre spans $(length(fibre)) dimensions")
  return 0.5 * dot(ξ, fibre.Λ * ξ) - dot(getInfovector(fibre), ξ) - fibre.logmass
end

energy(bundle::DensityBundlePoint, ξ::AbstractVector) = energy(bundle.fibre, ξ)

""" $SIGNATURES
Normalized conditional energy `−log p(ξ_F|ΔS)` with mean `m = v − W·ΔS`:
`φ = ½(ξ_F−m)ᵀΛ(ξ_F−m) − ½log|Λ| + (r/2)·log 2π`.
"""
function energy(c::GaussianConditional, ξ_F::AbstractVector, ΔS::AbstractVector)
  m = backsubstitute(c, ΔS)
  logdet_FF, r = _logDetPSD(c.Λ)
  return 0.5 * dot(ξ_F - m, c.Λ * (ξ_F - m)) - 0.5 * logdet_FF + 0.5 * r * log(2π)
end
