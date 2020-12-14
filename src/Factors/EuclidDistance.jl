
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


getSample(s::EuclidDistance, N::Int=1) = (reshape(rand(s.Z,N),:,N), )

function (s::EuclidDistance)( res::AbstractArray{<:Real},
                              fmd::FactorMetadata,
                              idx::Int,
                              meas::Tuple,
                              X1::AbstractArray{<:Real,2},
                              X2::AbstractArray{<:Real,2}  )
  #
  res[1] = meas[1][1,idx] - norm(X2[:,idx] - X1[:,idx])
  res[1] ^= 2
  return res[1]
end

# parametric specific functor
function (s::EuclidDistance{<:ParametricTypes})(X1::AbstractArray{<:Real},
                                                X2::AbstractArray{<:Real};
                                                userdata::Union{Nothing,FactorMetadata}=nothing )
  #
  # can I change userdata to a keyword arg, DF, No will be resolved with consolidation
  #
  # FIXME, replace if with dispatch
  if isa(s.Z, Normal)
    meas = mean(s.Z)
    σ = std(s.Z)
    # res = similar(X2)
    res = meas - norm(X2 - X1)
    res *= res
    return res/(σ^2)

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