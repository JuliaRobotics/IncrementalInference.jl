
"""
$(TYPEDEF)

Default prior on all dimensions of a variable node in the factor graph.  `Prior` is
not recommended when non-Euclidean dimensions are used in variables.
"""
struct Prior{T <: SamplableBelief} <: AbstractPrior
  Z::T
end
Prior(::UniformScaling) = Prior(Normal())

getManifold(pr::Prior) = TranslationGroup(getDimension(pr.Z))
# getSample(cf::CalcFactor{<:Prior}) = rand(cf.factor.Z)

# basic default
(s::CalcFactor{<:Prior})(z, x1) = z .- x1

## packed types are still developed by hand.  Future versions would likely use a @packable macro to write Protobuf safe versions of factors

"""
$(TYPEDEF)

Serialization type for Prior.
"""
Base.@kwdef mutable struct PackedPrior <: AbstractPackedFactor
  Z::PackedSamplableBelief
end
function convert(::Type{PackedPrior}, d::Prior)
  return PackedPrior(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{Prior}, d::PackedPrior)
  return Prior(convert(SamplableBelief, d.Z))
end

function DFG.pack(d::Prior)
  return PackedPrior(DFG.packDistribution(d.Z))
end
function DFG.unpack(d::PackedPrior)
  return Prior(DFG.unpackDistribution(d.Z))
end

#
