# agnostic factors for internal use

"""
$(TYPEDEF)
"""
mutable struct GenericMarginal <: RelativeObservation # AbstractRelativeRoots
  Zij::Array{Float64, 1}
  Cov::Array{Float64, 1}
  W::Array{Float64, 1}
  GenericMarginal() = new()
  GenericMarginal(a, b, c) = new(a, b, c)
end

getManifold(::GenericMarginal) = LieGroups.TranslationGroup(1)

getSample(::CalcFactor{<:GenericMarginal}) = [0]

Base.@kwdef mutable struct PackedGenericMarginal <: AbstractPackedFactor
  Zij::Array{Float64, 1}
  Cov::Array{Float64, 1}
  W::Array{Float64, 1}
end
function convert(::Type{PackedGenericMarginal}, d::GenericMarginal)
  return PackedGenericMarginal(d.Zij, d.Cov, d.W)
end
function convert(::Type{GenericMarginal}, d::PackedGenericMarginal)
  return GenericMarginal(d.Zij, d.Cov, d.W)
end

# ------------------------------------------------------------
