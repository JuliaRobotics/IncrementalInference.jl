
export LinearRelative, PackedLinearRelative



"""
$(TYPEDEF)

Default linear offset between two scalar variables.
"""
struct LinearRelative{N, T <: SamplableBelief} <: AbstractRelativeFactor
  Z::T
end

function LinearRelative{N}(::UniformScaling) where N
  newval = MvNormal(zeros(N), diagm(ones(N)))
  LinearRelative{N,typeof(newval)}(newval)
end
# LinearRelative(n::Int=1) = LinearRelative{n}()
LinearRelative(iden::UniformScaling=LinearAlgebra.I) = LinearRelative{1}(iden)
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
                                #can I change userdata to a keyword arg
  #
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
Serialization type for `LinearConditional` binary factor.
"""
mutable struct PackedLinearRelative <: PackedInferenceType
  Z::String
  PackedLinearRelative() = new()
  PackedLinearRelative(z::AS) where {AS <: AbstractString} = new(z)
end
function convert(::Type{PackedLinearRelative}, d::LinearRelative)
  PackedLinearRelative(string(d.Z))
end
function convert(::Type{LinearRelative}, d::PackedLinearRelative)
  LinearRelative(extractdistribution(d.Z))
end


# """
# $(TYPEDEF)
# Serialization type for `MixtureLinearConditional`.
# """
# mutable struct PackedMixtureLinearConditional <: PackedInferenceType
#   strs::Vector{String}
#   cat::String
#   PackedMixtureLinearConditional() = new()
#   PackedMixtureLinearConditional(z::Vector{<:AbstractString}, cstr::AS) where {AS <: AbstractString} = new(z, cstr)
# end
# function convert(::Type{PackedMixtureLinearConditional}, d::MixtureLinearConditional)
#   PackedMixtureLinearConditional(string.(d.Z), string(d.C))
# end
# function convert(::Type{MixtureLinearConditional}, d::PackedMixtureLinearConditional)
#   MixtureLinearConditional(extractdistribution.(d.strs), extractdistribution(d.cat))
# end



#