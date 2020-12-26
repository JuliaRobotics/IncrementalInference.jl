

"""
$(TYPEDEF)

Partial prior belief (absolute data) on any variable, given `<:SamplableBelief` and which dimensions of the intended variable.
"""
struct PartialPrior{T <: SamplableBelief,P <: Tuple} <: AbstractPrior
  Z::T
  partial::P
end
getSample(s::PartialPrior, N::Int=1) = (reshape(rand(s.Z,N),:,N), )



"""
$(TYPEDEF)

Serialization type for `PartialPrior`.
"""
mutable struct PackedPartialPrior <: PackedInferenceType
  Z::String
  partials::Vector{Int}
end

function convert(::Type{PackedPartialPrior}, d::PartialPrior)
  PackedPartialPrior(convert(PackedSamplableBelief, d.Z), [d.partial...;])
end
function convert(::Type{PartialPrior}, d::PackedPartialPrior)
  PartialPrior(convert(SamplableBelief, d.Z),(d.partials...,))
end



