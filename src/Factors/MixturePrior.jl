
export MixturePrior, PackedMixturePrior


_defaultNamesMixtures(N::Int) = ((Symbol[Symbol("c$i") for i in 1:N])...,)

"""
$(TYPEDEF)

Define a categorical mixture of prior beliefs on a variable (updated version).
"""
struct MixturePrior{N, S, T} <: AbstractPrior
  components::NamedTuple{S, T}
  diversity::Distributions.Categorical
  #derived values
  dims::Int
  labels::Vector{Int}
  # smpls::Array{Float64,2}
end

"""
    $SIGNATURES

Construction helper functions for any `MixturePrior`.
"""
MixturePrior( z::NamedTuple{S,T}, 
              c::Distributions.DiscreteNonParametric,
              N::Int=length(z) ) where {S, T} = MixturePrior{N,S,T}(z, c, size( rand(z[1],1), 1), zeros(Int, 0))
#
MixturePrior(z::NamedTuple, c::AbstractVector{<:Real} ) = MixturePrior(z, Distributions.Categorical(c))
MixturePrior(z::NamedTuple, c::NTuple{N, <:Real} ) where N = MixturePrior(z, [c...])
MixturePrior(z::NTuple{N,<:SamplableBelief}, c::Union{<:Distributions.DiscreteNonParametric, NTuple{N,<:Real}, <:AbstractVector{<:Real}} ) where N = MixturePrior(NamedTuple{_defaultNamesMixtures(N)}((z...,)), c)
MixturePrior(z::AbstractVector{<:SamplableBelief},c::Union{<:Distributions.DiscreteNonParametric, <:AbstractVector{<:Real}, NTuple{N,<:Real}} ) where N = MixturePrior((z...,), c)


# convert(PackedInferenceType, MixedPrior())

"""
$(TYPEDEF)

Serialization type for `MixedPrior`.
"""
mutable struct PackedMixturePrior <: PackedInferenceType
  components::Vector{String}
  diversity::String
  PackedMixturePrior() = new()
  PackedMixturePrior(z::Vector{<:AbstractString}, cstr::AS) where {AS <: AbstractString} = new(z, cstr)
end
function convert(::Type{PackedMixturePrior}, d::MixturePrior)

  PackedMixturePrior(string.(d.Z), string(d.C))
end
function convert(::Type{MixturePrior}, d::PackedMixturePrior)
  MixturePrior(extractdistribution.(d.components), extractdistribution(d.diversity))
end





#