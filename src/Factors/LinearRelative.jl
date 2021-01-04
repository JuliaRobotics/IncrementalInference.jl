
export LinearRelative, PackedLinearRelative



"""
$(TYPEDEF)

Default linear offset between two scalar variables.

```math
X_2 = X_1 + η_Z
```
"""
struct LinearRelative{N, T <: SamplableBelief} <: AbstractRelativeRoots
  Z::T
end


# need several helper constructors since the dimension over which LinearRelative will be used is unknown at this point
function LinearRelative{N}( z0::T=MvNormal(zeros(N), diagm(ones(N))) ) where {N, T <: SamplableBelief}
  #
  LinearRelative{N, T}(z0)
end

LinearRelative(::UniformScaling=LinearAlgebra.I) = LinearRelative{1}(MvNormal(zeros(1), diagm(ones(1))))
LinearRelative(nm::Distributions.ContinuousUnivariateDistribution) = LinearRelative{1, typeof(nm)}(nm)
LinearRelative(nm::MvNormal) = LinearRelative{length(nm.μ), typeof(nm)}(nm)
LinearRelative(nm::BallTreeDensity) = LinearRelative{Ndim(nm), typeof(nm)}(nm)

getDimension(::InstanceType{LinearRelative{N,<:SamplableBelief}}) where {N} = N
getManifolds(::InstanceType{LinearRelative{N,<:SamplableBelief}}) where {N} = tuple([:Euclid for i in 1:N]...)

getDomain(::InstanceType{LinearRelative{N,<:SamplableBelief}}) where N = ContinuousEuclid{N}
# getManifolds(fctType::Type{LinearRelative}) = getManifolds(getDomain(fctType))


getSample(cf::CalcFactor{<:LinearRelative}, N::Int=1) = (reshape(rand(cf.factor.Z,N),:,N), )



# new and simplified interface for both nonparametric and parametric
function (s::CalcFactor{<:LinearRelative})( res::AbstractVector{<:Real},
                                                  z,
                                                  x1,
                                                  x2  ) # where {M<:FactorMetadata,P<:Tuple,X<:AbstractVector}
  #
  # TODO convert to distance(distance(x2,x1),z) # or use dispatch on `-` -- what to do about `.-`
  res .= z - (x2 - x1)
  nothing
end




"""
$(TYPEDEF)
Serialization type for `LinearRelative` binary factor.
"""
mutable struct PackedLinearRelative <: PackedInferenceType
  Z::String
  PackedLinearRelative() = new()
  PackedLinearRelative(z::AS) where {AS <: AbstractString} = new(z)
end
function convert(::Type{PackedLinearRelative}, d::LinearRelative)
  PackedLinearRelative(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{LinearRelative}, d::PackedLinearRelative)
  LinearRelative(convert(SamplableBelief, d.Z))
end




#