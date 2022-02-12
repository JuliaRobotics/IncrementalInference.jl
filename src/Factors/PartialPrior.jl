

"""
$(TYPEDEF)

Partial prior belief (absolute data) on any variable, given `<:SamplableBelief` and which dimensions of the intended variable.

Notes
- If using [`AMP.ManifoldKernelDensity`](@ref), don't double partial.  Only define the partial in this `PartialPrior` container. 
  - Future TBD, consider using `AMP.getManifoldPartial` for more general abstraction.
"""
struct PartialPrior{T <: SamplableBelief,P <: Tuple} <: AbstractPrior
  Z::T
  partial::P
end

# TODO, standardize, but shows error on testPartialNH.jl
getSample(cf::CalcFactor{<:PartialPrior}) = samplePoint(cf.factor.Z)   # remove in favor of ManifoldSampling.jl
# getManifold(pp::PartialPrior) = TranslationGroup(length(pp.partial)) # uncomment


getManifold(pp::PartialPrior{<:PackedManifoldKernelDensity}) = pp.Z.manifold




"""
$(TYPEDEF)

Serialization type for `PartialPrior`.
"""
Base.@kwdef struct PackedPartialPrior <: AbstractPackedFactor
  Z::PackedSamplableBelief
  partials::Vector{Int}
end

function convert(::Type{PackedPartialPrior}, d::PartialPrior)
  PackedPartialPrior(convert(PackedSamplableBelief, d.Z), [d.partial...;])
end
function convert(::Type{PartialPrior}, d::PackedPartialPrior)
  PartialPrior(convert(SamplableBelief, d.Z),(d.partials...,))
end



