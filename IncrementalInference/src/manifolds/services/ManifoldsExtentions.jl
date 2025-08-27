
## ================================================================================================
## Manifold and ManifoldDiff use with Optim
## ================================================================================================

# Modified from: https://gist.github.com/mateuszbaran/0354c0edfb9cdf25e084a2b915816a09
"""
    ManifoldWrapper{TM<:AbstractManifold} <: Optim.Manifold
    
Adapts Manifolds.jl manifolds for use in Optim.jl
"""
struct ManifoldWrapper{TM <: AbstractManifold} <: Optim.Manifold
  M::TM
end

function Optim.retract!(M::ManifoldWrapper, x)
  ManifoldsBase.embed_project!(M.M, x, x)
  return x
end

function Optim.project_tangent!(M::ManifoldWrapper, g, x)
  ManifoldsBase.embed_project!(M.M, g, x, g)
  return g
end

## ================================================================================================
## AbstractPowerManifold with N as field to avoid excessive compiling time.
## ================================================================================================
struct NPowerManifold{𝔽, M <: AbstractManifold{𝔽}} <:
       AbstractPowerManifold{𝔽, M, NestedReplacingPowerRepresentation}
  manifold::M
  N::Int
end

Manifolds.get_iterator(M::NPowerManifold) = Base.OneTo(M.N)

function Manifolds.manifold_dimension(M::NPowerManifold)
  return manifold_dimension(M.manifold) * M.N
end

function Manifolds.get_vector!(M::NPowerManifold, Y, p, c, B::AbstractBasis)
  dim = manifold_dimension(M.manifold)
  rep_size = Manifolds.representation_size(M.manifold)
  v_iter = 1
  for i in Manifolds.get_iterator(M)
    Y[i] = get_vector(
      M.manifold,
      Manifolds._read(M, rep_size, p, i),
      # view(c, v_iter:(v_iter + dim - 1)),
      SVector{dim}(view(c, v_iter:(v_iter + dim - 1))),
      B,
    )
    v_iter += dim
  end
  return Y
end

function Manifolds.exp!(M::NPowerManifold, q, p, X)
  rep_size = Manifolds.representation_size(M.manifold)
  for i in Manifolds.get_iterator(M)
    q[i] = exp(
      M.manifold,
      Manifolds._read(M, rep_size, p, i),
      Manifolds._read(M, rep_size, X, i),
    )
  end
  return q
end

function Manifolds.compose!(M::NPowerManifold, x, p, q)
  rep_size = representation_size(M.manifold)
  for i in Manifolds.get_iterator(M)
    x[i] = compose(
      M.manifold,
      Manifolds._read(M, rep_size, p, i),
      Manifolds._read(M, rep_size, q, i),
    )
  end
  return x
end

function Manifolds.allocate_result(M::NPowerManifold, f, x...)
  if length(x) == 0
    return [Manifolds.allocate_result(M.manifold, f) for _ in Manifolds.get_iterator(M)]
  else
    return copy(x[1])
  end
end

function Manifolds.allocate_result(::NPowerManifold, ::typeof(get_vector), p, X)
  return copy(p)
end

## ================================================================================================
## ArrayPartition getPointIdentity (identity_element)
## ================================================================================================
# NOTE This will be removed once moved upstream to Manifolds.jl

import DistributedFactorGraphs: getPointIdentity

# fallback 
function DFG.getPointIdentity(G::AbstractLieGroup, ::Type{T} = Float64) where {T <: Real}
  # error("getPointIdentity not implemented for $G.")
  @warn("getPointIdentity not implemented on $G, falling back to identity_element")
  return identity_element(G)
end

function DFG.getPointIdentity(M::ProductManifold, ::Type{T} = Float64) where {T <: Real}
  return ArrayPartition(map(x -> getPointIdentity(x, T), M.manifolds))
end

function DFG.getPointIdentity(G::ProductGroup, ::Type{T} = Float64) where {T <: Real}
  M = G.manifold
  return ArrayPartition(map(x -> getPointIdentity(x, T), M.manifolds))
end

function DFG.getPointIdentity(
  G::LieGroups.TranslationGroup{ℝ, TypeParameter{Tuple{N}}},
  ::Type{T} = Float64,
) where {N, T <: Real}
  return zeros(SVector{N, T})
end

#TODO test
function DFG.getPointIdentity(
  PrG::LieGroup{𝔽, Op, M},
  ::Type{T} = Float64,
) where {𝔽, Op <: AbstractProductGroupOperation, M <: ProductManifold, T <: Real}
  PrM = PrG.manifold
  ε = map(G -> getPointIdentity(G, T), map(LieGroup, PrM.manifolds, PrG.op.operations))
  return ArrayPartition(ε)
end

function DFG.getPointIdentity(
  VG::LieGroups.ValidationLieGroup,
  ::Type{T} = Float64,
) where {T <: Real}
  G = VG.lie_group
  return LieGroups.ValidationMPoint(getPointIdentity(G, T))
end

function DFG.getPointIdentity(
  @nospecialize(M::PowerManifold),
  ::Type{T} = Float64,
) where {T <: Real}
  N = Manifolds.get_iterator(M).stop
  return fill(getPointIdentity(M.manifold, T), N)
end

function DFG.getPointIdentity(M::NPowerManifold, ::Type{T} = Float64) where {T <: Real}
  return fill(getPointIdentity(M.manifold, T), M.N)
end

function DFG.getPointIdentity(
  ::typeof(SpecialEuclideanGroup(2; variant = :right)),
  ::Type{T} = Float64,
) where {T <: Real}
  N = 2
  return ArrayPartition(zeros(SVector{N, T}), SMatrix{N, N, T}(I))
end

function DFG.getPointIdentity(
  ::typeof(SpecialEuclideanGroup(3; variant = :right)),
  ::Type{T} = Float64,
) where {T <: Real}
  N = 3
  return ArrayPartition(zeros(SVector{N, T}), SMatrix{N, N, T}(I))
end

function DFG.getPointIdentity(
  G::SpecialOrthogonalGroup{TypeParameter{Tuple{N}}},
  ::Type{T} = Float64,
) where {N, T <: Real}
  return SMatrix{N, N, T}(I)
end

function DFG.getPointIdentity(
  G::LieGroups.TranslationGroup{TypeParameter{Tuple{N}}},
  ::Type{T} = Float64,
) where {N, T <: Real}
  return zeros(SVector{N, T})
end

function DFG.getPointIdentity(G::AbstractLieGroup{ℝ,AdditionGroupOperation,<:Circle{ℝ}}, ::Type{T} = Float64) where {T <: Real}
  return [zero(T)] #FIXME we cannot support scalars yet
end

function DFG.getPointIdentity(G::typeof(LieGroups.CircleGroup()), ::Type{T} = Float64) where {T <: Real}
  return fill(Complex{T}(1,0))
end
