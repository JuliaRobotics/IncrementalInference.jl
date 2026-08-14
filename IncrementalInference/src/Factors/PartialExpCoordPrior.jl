## ======================================================================================
## Partial prior on a subspace of the Lie algebra (tangent space at the identity).
##
## Preferred over `PartialPrior` for Lie group variables: the observed subspace is stated explicitly
## in the Lie algebra rather than implied by index positions on whatever representation the variable
## happens to use.
## ======================================================================================
# example of 2 `e1 = [1, 0, 0]` "partials" shows why:
# - `f(R) = vee(log(R))ᵀ * e₁  # SO(3) -> ℝ¹`
# - `f(R) = R * e₁             # SO(3) -> S²`

"""
    $TYPEDEF

A partial prior constraining a subspace of the Lie algebra of a group variable.

The variable is mapped into the Lie algebra (the tangent space at the identity) and expressed in
coordinates by `vee` (`DefaultLieAlgebraOrthogonalBasis`); the observed subspace is then selected by `U`:

    r = z - Uᵀ · vee(LieAlgebra(G), log(G, g))

`U` is `n × k` with orthonormal columns spanning the observed subspace, `n = manifold_dimension(G)`.
`getManifold` returns the `k`-dimensional *residual* space, so the factor's dimension is the number of
observed directions rather than the variable's.

# Constructors

    PartialExpCoordPrior(G, Z, U::AbstractMatrix)   # arbitrary observed subspace
    PartialExpCoordPrior(G, Z, partial::Tuple)      # observe individual coordinates

The tuple form builds the corresponding selection matrix, so observing coordinates is just the
axis-aligned special case of observing a subspace.

`log` is taken **about the group identity**, so the coordinates are those of the Lie algebra in the
basis `vee` uses. The factor is a locally valid submersion on any Lie group provided the variable
stays inside the injectivity radius, where that parameterization is well defined.

$(TYPEDFIELDS)
"""
struct PartialExpCoordPrior{
  G <: LieGroups.AbstractLieGroup,
  T <: SamplableBelief,
  UT <: AbstractMatrix,
} <: AbstractPriorObservation
  G::G
  Z::T
  U::UT
end

"""
    $SIGNATURES

Observe individual Lie algebra coordinates — the axis-aligned special case, which builds the
corresponding selection matrix for `U`.
"""
function PartialExpCoordPrior(
  G::LieGroups.AbstractLieGroup,
  Z::SamplableBelief,
  partial::Tuple,
)
  n = manifold_dimension(G)
  U = zeros(n, length(partial))
  for (j, i) in enumerate(partial)
    @assert 1 <= i <= n "partial index $i outside 1:$n for $G"
    U[i, j] = 1.0
  end
  return PartialExpCoordPrior(G, Z, U)
end

DFG.getManifold(prior::PartialExpCoordPrior) = LieGroups.TranslationGroup(size(prior.U, 2))

function (cf::CalcFactor{<:PartialExpCoordPrior})(z, x1)
  G = cf.factor.G
  Xc = vee(LieAlgebra(G), log(G, x1))
  return z .- cf.factor.U' * Xc
end
