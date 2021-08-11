
export EuclidDistance, PackedEuclidDistance



"""
$(TYPEDEF)

Default linear offset between two scalar variables.
"""
struct EuclidDistance{T <: SamplableBelief} <: AbstractRelativeMinimize
  Z::T
end

EuclidDistance(::UniformScaling=LinearAlgebra.I) = EuclidDistance(Normal())

# consider a different group?
getManifold(::InstanceType{EuclidDistance}) = TranslationGroup(1)
getDimension(::InstanceType{<:EuclidDistance}) = 1


function getSample(cf::CalcFactor{<:EuclidDistance}, N::Int=1)
  ret = [rand(cf.factor.Z,1) for _ in 1:N]
  (ret, )
end

# new and simplified interface for both nonparametric and parametric
function (s::CalcFactor{<:EuclidDistance})(z, x1, x2)
  # @info "distance?" z x1 x2
  return z .- norm(x2 .- x1)
end


Base.convert(::Type{<:MB.AbstractManifold}, ::InstanceType{EuclidDistance}) = Manifolds.TranslationGroup(1)


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