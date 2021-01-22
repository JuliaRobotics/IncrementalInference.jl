
export EuclidDistance, PackedEuclidDistance



"""
$(TYPEDEF)

Default linear offset between two scalar variables.
"""
struct EuclidDistance{T <: SamplableBelief} <: AbstractRelativeMinimize
  Z::T
end

EuclidDistance(::UniformScaling=LinearAlgebra.I) = EuclidDistance(Normal())

getDimension(::InstanceType{EuclidDistance{<:SamplableBelief}}) = 1

# before consolidation, see RoME.jl #244
getManifolds(::InstanceType{EuclidDistance{<:SamplableBelief}}) = (:Euclid,)
getDomain(::InstanceType{EuclidDistance{<:SamplableBelief}}) = ContinuousEuclid{1}


getSample(cf::CalcFactor{<:EuclidDistance}, N::Int=1) = (reshape(rand(cf.factor.Z,N),:,N), )


# new and simplified interface for both nonparametric and parametric
function (s::CalcFactor{<:EuclidDistance})( residual::AbstractVector{<:Real},
                                            z,
                                            x1,
                                            x2 )
  #
  residual .= z .- norm(x2 .- x1)

  # v0.21+, should not return cost, but rather populate residual
  nothing
end




"""
$(TYPEDEF)
Serialization type for `EuclidDistance` binary factor.
"""
mutable struct PackedEuclidDistance <: PackedInferenceType
  mimeType::String
  Z::String
  # PackedEuclidDistance() = new()
  # PackedEuclidDistance(z::AS) where {AS <: AbstractString} = new(z)
end

function convert(::Type{PackedEuclidDistance}, d::EuclidDistance)
  PackedEuclidDistance( "/application/JuliaLang/PackedSamplableBelief",
                        convert(PackedSamplableBelief, d.Z) )
end

function convert(::Type{<:EuclidDistance}, d::PackedEuclidDistance)
  EuclidDistance(convert(SamplableBelief, d.Z))
end







#