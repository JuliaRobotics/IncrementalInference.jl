
export EuclidDistance, PackedEuclidDistance

"""
$(TYPEDEF)

Default linear offset between two scalar variables.
"""
struct EuclidDistance{T <: SamplableBelief} <: AbstractManifoldMinimize # AbstractRelativeMinimize
  Z::T
end

EuclidDistance(::UniformScaling = LinearAlgebra.I) = EuclidDistance(Normal())

# consider a different group?
getManifold(::InstanceType{EuclidDistance}) = TranslationGroup(1)
getDimension(::InstanceType{<:EuclidDistance}) = 1

# new and simplified interface for both nonparametric and parametric
(s::CalcFactor{<:EuclidDistance})(z, x1, x2) = z .- norm(x2 .- x1)

function Base.convert(::Type{<:MB.AbstractManifold}, ::InstanceType{EuclidDistance})
  return Manifolds.TranslationGroup(1)
end

"""
$(TYPEDEF)
Serialization type for `EuclidDistance` binary factor.
"""
Base.@kwdef mutable struct PackedEuclidDistance <: AbstractPackedFactor
  _type::String
  Z::PackedSamplableBelief
end

function convert(::Type{PackedEuclidDistance}, d::EuclidDistance)
  return PackedEuclidDistance(
    "/application/JuliaLang/PackedSamplableBelief",
    convert(PackedSamplableBelief, d.Z),
  )
end

function convert(::Type{<:EuclidDistance}, d::PackedEuclidDistance)
  return EuclidDistance(convert(SamplableBelief, d.Z))
end

#
