
export EuclidDistance, PackedEuclidDistance



"""
$(TYPEDEF)

Default linear offset between two scalar variables.
"""
struct EuclidDistance{T <: SamplableBelief} <: AbstractRelativeMinimize
  Z::T
end

EuclidDistance(::UniformScaling=LinearAlgebra.I) = EuclidDistance(Normal())

getDimension(::InstanceType{<:EuclidDistance}) = 1

# before consolidation, see RoME.jl #244
getManifolds(::InstanceType{<:EuclidDistance}) = (:Euclid,)
getDomain(::InstanceType{<:EuclidDistance}) = ContinuousEuclid{1}


getSample(cf::CalcFactor{<:EuclidDistance}, N::Int=1) = (reshape(rand(cf.factor.Z,N),1,N), )


# new and simplified interface for both nonparametric and parametric
function (s::CalcFactor{<:EuclidDistance})(z, x1, x2)
  # v0.21+, should return residual and not have residual parameter
  return z .- norm(x2 .- x1)
end


Base.convert(::Type{<:ManifoldsBase.Manifold}, ::InstanceType{EuclidDistance}) = Manifolds.Euclidean(1)


"""
$(TYPEDEF)
Serialization type for `EuclidDistance` binary factor.
"""
mutable struct PackedEuclidDistance <: PackedInferenceType
  _type::String
  Z::String
end

function convert(::Type{PackedEuclidDistance}, d::EuclidDistance)
  PackedEuclidDistance( "/application/JuliaLang/PackedSamplableBelief",
                        convert(PackedSamplableBelief, d.Z) )
end

function convert(::Type{<:EuclidDistance}, d::PackedEuclidDistance)
  EuclidDistance(convert(SamplableBelief, d.Z))
end







#