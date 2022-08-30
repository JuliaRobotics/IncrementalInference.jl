
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
function LinearRelative{N}(
  z0::T = MvNormal(zeros(N), diagm(ones(N))),
) where {N, T <: SamplableBelief}
  #
  return LinearRelative{N, T}(z0)
end

function LinearRelative(::UniformScaling = LinearAlgebra.I)
  return LinearRelative{1}(MvNormal(zeros(1), diagm(ones(1))))
end
function LinearRelative(nm::Distributions.ContinuousUnivariateDistribution)
  return LinearRelative{1, typeof(nm)}(nm)
end
LinearRelative(nm::MvNormal) = LinearRelative{length(nm.μ), typeof(nm)}(nm)
function LinearRelative(nm::Union{<:BallTreeDensity, <:ManifoldKernelDensity})
  return LinearRelative{Ndim(nm), typeof(nm)}(nm)
end

getManifold(::InstanceType{LinearRelative{N}}) where {N} = getManifold(ContinuousEuclid{N})

# TODO standardize
getDimension(::InstanceType{LinearRelative{N}}) where {N} = N

# new and simplified interface for both nonparametric and parametric
function (s::CalcFactor{<:LinearRelative})(z, x1, x2)
  # TODO convert to distance(distance(x2,x1),z) # or use dispatch on `-` -- what to do about `.-`
  return z .- (x2 .- x1)
end

function Base.convert(
  ::Type{<:MB.AbstractManifold},
  ::InstanceType{LinearRelative{N}},
) where {N}
  return Manifolds.TranslationGroup(N)
end

"""
$(TYPEDEF)
Serialization type for `LinearRelative` binary factor.
"""
Base.@kwdef mutable struct PackedLinearRelative <: AbstractPackedFactor
  Z::PackedSamplableBelief
end
function convert(::Type{PackedLinearRelative}, d::LinearRelative)
  return PackedLinearRelative(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{LinearRelative}, d::PackedLinearRelative)
  return LinearRelative(convert(SamplableBelief, d.Z))
end

#
