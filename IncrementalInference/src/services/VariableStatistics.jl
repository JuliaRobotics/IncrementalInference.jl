
function Statistics.mean(vartype::StateType, args...; kwargs...)
  return mean(getManifold(vartype), args...; kwargs...)
end
function Statistics.std(vartype::StateType, args...; kwargs...)
  return std(getManifold(vartype), args...; kwargs...)
end
function Statistics.var(vartype::StateType, args...; kwargs...)
  return var(getManifold(vartype), args...; kwargs...)
end

function Statistics.cov(
  vartype::StateType,
  ptsArr::AbstractVector;
  basis::Manifolds.AbstractBasis = Manifolds.DefaultOrthogonalBasis(),
  kwargs...,
)
  return cov(getManifold(vartype), ptsArr; basis, kwargs...)
end

#TODO check performance and FIXME on makemutalbe might not be needed any more
function calcStdBasicSpread(vartype::StateType, ptsArr::AbstractVector) # {P}) where {P}
  # _makemutable(s) = s
  # _makemutable(s::StaticArray{Tuple{S},T,N}) where {S,T,N} = MArray{Tuple{S},T,N,S}(s)
  # _makemutable(s::SMatrix{N,N,T,D}) where {N,T,D} = MMatrix{N,N,T,D}(s)
  
  # FIXME, silly conversion since Manifolds.std internally replicates eltype ptsArr which doesn't work on StaticArrays
  # σ = std(vartype, _makemutable.(ptsArr))

  μ = mean(vartype, ptsArr, GeodesicInterpolation())
  σ = std(vartype, ptsArr, μ)

  #if no std yet, set to 1
  msst = 1e-10 < σ ? σ : 1.0
  return msst
end

#TODO consolidate
function calcMeanCovar(vari::VariableCompute, solvekey = :default)
  pts = getState(vari, solvekey).val
  μ = mean(getManifold(vari), pts)
  Σ = cov(getStateKind(vari), pts)
  return μ, Σ
end
