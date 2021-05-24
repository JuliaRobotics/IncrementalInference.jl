
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
LinearRelative(nm::Union{<:BallTreeDensity,<:ManifoldKernelDensity}) = LinearRelative{Ndim(nm), typeof(nm)}(nm)

# getManifold(::InstanceType{LinearRelative{N}}) where {N} = Euclidean(N)
# getManifolds(::T) where {T <: LinearRelative} = convert(Tuple, getManifold(T))
# getManifolds(::Type{<:T}) where {T <: LinearRelative} = convert(Tuple, getManifold(T))
# getManifolds(fctType::Type{LinearRelative}) = getManifolds(getDomain(fctType))

getManifold(::InstanceType{LinearRelative{N}}) where N = ContinuousEuclid{N}

# TODO standardize
getDimension(::InstanceType{LinearRelative{N}}) where {N} = N

getSample(cf::CalcFactor{<:LinearRelative}, N::Int=1) = (reshape(rand(cf.factor.Z,N),:,N), )



# new and simplified interface for both nonparametric and parametric
function (s::CalcFactor{<:LinearRelative})(z, x1, x2) 
  # TODO convert to distance(distance(x2,x1),z) # or use dispatch on `-` -- what to do about `.-`
  # v0.21+, should return residual
  return z .- (x2 .- x1)
end



Base.convert(::Type{<:MB.AbstractManifold}, ::InstanceType{LinearRelative{N}}) where N = Manifolds.Euclidean(N) 



"""
$(TYPEDEF)
Serialization type for `LinearRelative` binary factor.
"""
mutable struct PackedLinearRelative <: PackedInferenceType
  Z::String
  # PackedLinearRelative() = new()
  # PackedLinearRelative(z::AS) where {AS <: AbstractString} = new(z)
end
function convert(::Type{PackedLinearRelative}, d::LinearRelative)
  PackedLinearRelative(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{LinearRelative}, d::PackedLinearRelative)
  LinearRelative(convert(SamplableBelief, d.Z))
end




#