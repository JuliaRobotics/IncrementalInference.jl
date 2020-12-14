
export LinearRelative, PackedLinearRelative



"""
$(TYPEDEF)

Default linear offset between two scalar variables.
"""
struct LinearRelative{N, T <: SamplableBelief} <: AbstractRelativeRoots
  Z::T

  # # Julia 1.5.2 still requires this inner constructor given all default helpers below?
  # LinearRelative{N, T}(z0::T) where {N, T <: SamplableBelief} = new{N,T}(z0)
end

# function LinearRelative{N, T}(::UniformScaling=LinearAlgebra.I,
#                               z0::T=MvNormal(zeros(N), diagm(ones(N))) ) where {N, T <: SamplableBelief}
#   #
#   LinearRelative{N,T}(z0)
# end

function LinearRelative{N}( z0::T=MvNormal(zeros(N), diagm(ones(N))) ) where {N, T <: SamplableBelief}
  #
  LinearRelative{N, T}(z0)
end

LinearRelative(::UniformScaling=LinearAlgebra.I) = LinearRelative{1}(MvNormal(zeros(1), diagm(ones(1))))
LinearRelative(nm::Distributions.ContinuousUnivariateDistribution) = LinearRelative{1, typeof(nm)}(nm)
LinearRelative(nm::MvNormal) = LinearRelative{length(nm.μ), typeof(nm)}(nm)
LinearRelative(nm::BallTreeDensity) = LinearRelative{Ndim(nm), typeof(nm)}(nm)

getDimension(::Type{LinearRelative{N,<:SamplableBelief}}) where {N} = N
getManifolds(::Type{LinearRelative{N,<:SamplableBelief}}) where {N} = tuple([:Euclid for i in 1:N]...)

getDomain(::InstanceType{LinearRelative}) = ContinuousScalar
getManifolds(fctType::Type{LinearRelative}) = getManifolds(getDomain(fctType))


getSample(s::LinearRelative, N::Int=1) = (reshape(rand(s.Z,N),:,N), )

function (s::LinearRelative)( res::AbstractArray{<:Real},
                              userdata::FactorMetadata,
                              idx::Int,
                              meas::Tuple,
                              X1::AbstractArray{<:Real,2},
                              X2::AbstractArray{<:Real,2}  )
  #
  res[:] = meas[1][:,idx] - (X2[:,idx] - X1[:,idx])
  nothing
end

# parametric specific functor
function (s::LinearRelative{N,<:ParametricTypes})(
                                X1::AbstractArray{<:Real},
                                X2::AbstractArray{<:Real};
                                userdata::Union{Nothing,FactorMetadata}=nothing ) where N
  #
  # can I change userdata to a keyword arg, DF, No will be resolved with consolidation
  # FIXME, replace if with dispatch
  if isa(s.Z, Normal)
    meas = mean(s.Z)
    σ = std(s.Z)
    # res = similar(X2)
    res = meas - (X2[1] - X1[1])
    return (res/σ) .^ 2

  elseif isa(s.Z, MvNormal)
    meas = mean(s.Z)
    iΣ = invcov(s.Z)
    #TODO confirm math : Σ^(1/2)*X
    res = meas .- (X2 .- X1)
    return res' * iΣ * res

  else
    #this should not happen
    @error("$s not suported, please use non-parametric")
  end
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