## ================================================================================================
## AbstractPowerManifold with N as field to avoid excessive compiling time.
## ================================================================================================
struct NPowerManifold{𝔽, M<:AbstractManifold{𝔽}} <: AbstractPowerManifold{𝔽, M, NestedReplacingPowerRepresentation}
  manifold::M
  N::Int
end

Manifolds.get_iterator(M::NPowerManifold) = Base.OneTo(M.N)

function Manifolds.manifold_dimension(M::NPowerManifold)
  return manifold_dimension(M.manifold) * M.N
end

function Manifolds.get_vector!(M::NPowerManifold, Y, p, c, B::AbstractBasis)
  dim = manifold_dimension(M.manifold)
  rep_size = representation_size(M.manifold)
  v_iter = 1
  for i in Manifolds.get_iterator(M)
      Y[i] = get_vector(
          M.manifold,
          Manifolds._read(M, rep_size, p, i),
          view(c,v_iter:(v_iter + dim - 1)),
          B,
      )
      v_iter += dim
  end
  return Y
end

function Manifolds.exp!(M::NPowerManifold, q, p, X)
  rep_size = representation_size(M.manifold)
  for i in Manifolds.get_iterator(M)
      q[i] = exp(M.manifold, Manifolds._read(M, rep_size, p, i), Manifolds._read(M, rep_size, X, i))
  end
  return q
end

## ================================================================================================
## ArrayPartition getPointIdentity (identity_element)
## ================================================================================================
# NOTE This will be removed once moved upstream to Manifolds.jl

import DistributedFactorGraphs: getPointIdentity


function getPointIdentity(G::ProductGroup, ::Type{T}=Float64) where T<:Real
  M = G.manifold
  return ArrayPartition(map(x->getPointIdentity(x, T), M.manifolds))
end

# fallback 
function getPointIdentity(G::GroupManifold,::Type{T}=Float64) where T<:Real
  return error("getPointIdentity not implemented on G")
end

function getPointIdentity(@nospecialize(G::ProductManifold),::Type{T}=Float64) where T<:Real
  return ArrayPartition(map(x->getPointIdentity(x,T), G.manifolds))
end

function getPointIdentity(@nospecialize(M::PowerManifold),::Type{T}=Float64) where T<:Real
  N = Manifolds.get_iterator(M).stop
  return fill(getPointIdentity(M.manifold, T), N)
end

function getPointIdentity(M::NPowerManifold,::Type{T}=Float64) where T<:Real
  return fill(getPointIdentity(M.manifold, T), M.N)
end

function getPointIdentity(G::SemidirectProductGroup,::Type{T}=Float64) where T<:Real
  M = base_manifold(G)
  N, H = M.manifolds
  np = getPointIdentity(N,T)
  hp = getPointIdentity(H,T)
  return ArrayPartition(np, hp)
end

#FIXME fix back to SA
function getPointIdentity(G::SpecialOrthogonal{N},::Type{T}=Float64) where {N,T<:Real}
  return SMatrix{N,N, T}(I)
  # return Matrix{T}(I, N, N)
end

function getPointIdentity(G::TranslationGroup{Tuple{N}},::Type{T}=Float64) where{N,T<:Real}
  return zeros(SVector{N,T})
  # return zeros(T,N)
end

function getPointIdentity(G::RealCircleGroup,::Type{T}=Float64) where T<:Real
    return zero(T)
    # return [zero(T)]
end
