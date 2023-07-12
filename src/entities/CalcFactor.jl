
abstract type AbstractMaxMixtureSolver end


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

Related 

[`CalcFactorMahalanobis`](@ref), [`CommonConvWrapper`](@ref)
"""
struct CalcFactor{
  FT <: AbstractFactor, 
  X, 
  C, 
  VT <: Tuple, 
  M <: AbstractManifold
}
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



"""
$TYPEDEF

Internal parametric extension to [`CalcFactor`](@ref) used for buffering measurement and calculating Mahalanobis distance

Related

[`CalcFactor`](@ref)
"""
struct CalcFactorMahalanobis{N, D, L, S <: Union{Nothing, AbstractMaxMixtureSolver}}
  faclbl::Symbol
  calcfactor!::CalcFactor
  varOrder::Vector{Symbol}
  meas::NTuple{N, <:AbstractArray}
  iΣ::NTuple{N, SMatrix{D, D, Float64, L}}
  specialAlg::S
end




struct CalcFactorManopt{
  D,
  L,
  FT <: AbstractFactor,
  M <: AbstractManifold,
  MEAS <: AbstractArray,
}
  faclbl::Symbol
  calcfactor!::CalcFactor{FT, Nothing, Nothing, Tuple{}, M}
  varOrder::Vector{Symbol}
  varOrderIdxs::Vector{Int}
  meas::MEAS
  iΣ::SMatrix{D, D, Float64, L}
  sqrt_iΣ::SMatrix{D, D, Float64, L}
end


