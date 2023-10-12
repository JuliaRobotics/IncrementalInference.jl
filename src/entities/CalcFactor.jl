
abstract type AbstractMaxMixtureSolver end


abstract type CalcFactor{T<:AbstractFactor} end


"""
$TYPEDEF

User factor interface method for computing the residual values of factors.

Notes
- Also see #467 on API consolidation

```julia
function (cf::CalcFactor{<:LinearRelative})(res::AbstractVector{<:Real}, z, xi, xj)
  cf.variablelist
  cf.cache
  # generic on-manifold residual function 
  return distance(z, distance(xj, xi))
end
```

DevNotes
- Follow the Github project in IIF to better consolidate CCW FMD CPT CF CFM
- TODO CalcFactorNormSq is a step towards having a dedicated structure for non-parametric solve.
  CalcFactorNormSq will calculate the Norm Squared of the factor.
Related 

[`CalcFactorMahalanobis`](@ref), [`CommonConvWrapper`](@ref)
"""
struct CalcFactorNormSq{
  FT <: AbstractFactor, 
  X, 
  C, 
  VT <: Tuple, 
  M <: AbstractManifold
} <: CalcFactor{FT}
  """ the interface compliant user object functor containing the data and logic """
  factor::FT
  """ what is the sample (particle) id for which the residual is being calculated """
  _sampleIdx::Int
  """ legacy support for variable values old functor residual functions.
      TBD, this is still being used by DERelative factors. """
  _legacyParams::X
  """ allow threading for either sampling or residual calculations (workaround for thread yield issue) """
  _allowThreads::Bool
  """ user cache of arbitrary type, overload the [`preambleCache`](@ref) function. NOT YET THREADSAFE """
  cache::C

  ## TODO Consolidation WIP with FactorMetadata
  # full list of variables connected to the factor
  # TODO make sure this list is of the active hypo only
  fullvariables::VT # Vector{<:DFGVariable} # FIXME change to tuple for better type stability
  # which index is being solved for?
  solvefor::Int
  manifold::M
end

#TODO deprecate after CalcFactor is updated to CalcFactorNormSq
function CalcFactor(args...; kwargs...)
  Base.depwarn(
    "`CalcFactor` changed to an abstract type, use CalcFactorNormSq, CalcFactorMahalanobis, or CalcFactorResidual",
    :CalcFactor
  )
  
  CalcFactorNormSq(args...; kwargs...)
end

"""
$TYPEDEF

Internal parametric extension to [`CalcFactor`](@ref) used for buffering measurement and calculating Mahalanobis distance

Related

[`CalcFactor`](@ref)
"""
struct CalcFactorMahalanobis{
  FT,
  N,
  C,
  MEAS<:AbstractArray,
  D,
  L,
  S <: Union{Nothing, AbstractMaxMixtureSolver}
} <: CalcFactor{FT}
  faclbl::Symbol
  factor::FT
  cache::C
  varOrder::Vector{Symbol}
  meas::NTuple{N, MEAS}
  iΣ::NTuple{N, SMatrix{D, D, Float64, L}}
  specialAlg::S
end


struct CalcFactorResidual{
  FT <: AbstractFactor,
  C,
  D,
  L,
  P,
  MEAS <: AbstractArray,
  N
} <: CalcFactor{FT}
  faclbl::Symbol
  factor::FT
  cache::C
  varOrder::NTuple{N, Symbol}
  varOrderIdxs::NTuple{N, Int}
  points::P #TODO remove or not?
  meas::MEAS
  iΣ::SMatrix{D, D, Float64, L} #TODO remove or not?
  sqrt_iΣ::SMatrix{D, D, Float64, L}
end

_nvars(::CalcFactorResidual{FT, C, D, L, P, MEAS, N}) where {FT, C, D, L, P, MEAS, N} = N
# _typeof_meas(::CalcFactorManopt{FT, C, D, L, MEAS, N}) where {FT, C, D, L, MEAS, N} = MEAS
DFG.getDimension(::CalcFactorResidual{FT, C, D, L, P, MEAS, N}) where {FT, C, D, L, P, MEAS, N} = D

# workaround for issue #1781
import Base: getproperty
function Base.getproperty(cf::CalcFactor, f::Symbol)
  if f === :manifold
    # assumes constant propagation to avoid allocations in residual functions getManifold(factor)
    getManifold(cf.factor)
  else
    getfield(cf, f)
  end
end
