

"""
$(TYPEDEF)

Define a categorical mixture of (relative) likelihood beliefs between any two variables.
"""
struct MixtureRelative{N, S, T<:Tuple} <: AbstractRelativeFactor
  components::NamedTuple{S,T}
  diversity::Distributions.Categorical
  dims::Int
  labels::Vector{Int}
end

# NamedTuple{_defaultNamesMixtures(N)}((z...,)), c

MixtureRelative(z::NamedTuple{S,T}, 
                c::Distributions.DiscreteNonParametric ) where {S, T} = MixtureRelative{length(z),S,T}(z, c, size( rand(z[1],1), 1), zeros(Int, 0))
#
MixtureRelative(z::NamedTuple{S,T}, c::AbstractVector{<:Real}) where {S,T} = MixtureRelative(z, Categorical([c...]) )
MixtureRelative(z::NamedTuple{S,T}, c::NTuple{N,<:Real}) where {N,S,T} = MixtureRelative(z, [c...] )
MixtureRelative(z::AbstractVector{<:SamplableBelief}, c::Union{<:Distributions.DiscreteNonParametric, <:AbstractVector{<:Real}, <:NTuple{N,<:Real}} ) where N = MixtureRelative(NamedTuple{_defaultNamesMixtures(length(z))}((z...,)), c)
MixtureRelative(z::Tuple, c::Union{<:Distributions.DiscreteNonParametric, <:AbstractVector{<:Real}, <:NTuple{N,<:Real}} ) where N = MixtureRelative(NamedTuple{_defaultNamesMixtures(length(z))}(z), c)


# getSample(s::MixtureRelative, N::Int=1) = (reshape.(rand.(s.Z, N),1,:)..., rand(s.C, N))


# function (s::MixtureLinearConditional)( res::AbstractArray{<:Real},
#                                         userdata::FactorMetadata,
#                                         idx::Int,
#                                         meas::Tuple,
#                                         X1::AbstractArray{<:Real,2},
#                                         X2::AbstractArray{<:Real,2}  )
#   #
#   res[1] = meas[meas[end][idx]][idx] - (X2[1,idx] - X1[1,idx])
#   nothing
# end



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







# struct MixtureLinearConditional{T} <: AbstractRelativeFactor



  #