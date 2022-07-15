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

function getPointIdentity(G::ProductManifold,::Type{T}=Float64) where T<:Real
  return ArrayPartition(map(x->getPointIdentity(x,T), G.manifolds))
end

function getPointIdentity(M::Manifolds.PowerManifoldNestedReplacing,::Type{T}=Float64) where T<:Real
  N = Manifolds.get_iterator(M).stop
  return fill(getPointIdentity(M.manifold, T), N)
end

function getPointIdentity(M::PowerManifold,::Type{T}=Float64) where T<:Real
  N = Manifolds.get_iterator(M).stop
  return fill(getPointIdentity(M.manifold, T), N)
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
  # return SMatrix{N,N, T}(I)
  return Matrix{T}(I, N, N)
end

function getPointIdentity(G::TranslationGroup{Tuple{N}},::Type{T}=Float64) where{N,T<:Real}
  # return zeros(SVector{N,T})
  return zeros(T,N)
end

function getPointIdentity(G::RealCircleGroup,::Type{T}=Float64) where T<:Real
    # return zero(T)
    return [zero(T)]
end
