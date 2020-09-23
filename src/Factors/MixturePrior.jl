


"""
$(TYPEDEF)

Define a categorical mixture of prior beliefs on a variable (updated version).
"""
struct MixturePrior{S, T <: Tuple} <: AbstractPrior
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
              c::Distributions.DiscreteNonParametric ) where {S, T} = MixturePrior{S,T}(z, c, size( rand(z[1],1), 1), zeros(Int, 0))
#
MixturePrior(z::NamedTuple, c::AbstractVector{<:Real} ) = MixturePrior(z, Distributions.Categorical(c))
MixturePrior(z::NamedTuple, c::NTuple{N, <:Real} ) where N = MixturePrior(z, [c...])
MixturePrior(z::NTuple{N,<:SamplableBelief}, c::Union{<:Distributions.DiscreteNonParametric, NTuple{N,<:Real}, <:AbstractVector{<:Real}} ) where N = MixturePrior(NamedTuple{((Symbol[Symbol("c$i") for i in 1:N])...,)}((z...,)), c)
MixturePrior(z::AbstractVector{<:SamplableBelief},c::Union{<:Distributions.DiscreteNonParametric, <:AbstractVector{<:Real}, NTuple{N,<:Real}} ) where N = MixturePrior((z...,), c)



function Base.resize!(mp::MixturePrior, s::Int)
  resize!(mp.labels, s)
end

# TODO make in-place memory version
function getSample(s::MixturePrior, N::Int=1)
  #out memory should be right size first
  (length(s.labels) != N) && resize!(s, N)
  s.labels .= rand(s.diversity, N)
  smpls = Array{Float64,2}(undef,s.dims,N)
  for i in 1:N
    mixComponent = s.components[s.labels[i]]
    smpls[:,i] .= rand(mixComponent,1)
  end
  (smpls, s.labels)
end



MixturePrior{T}(x...) where T = @error("`MixturePrior{$T}` is deprecated, use the new `NamedTuple`-based `MixturePrior{S,T}` instead.")


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