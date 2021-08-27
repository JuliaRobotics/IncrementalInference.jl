

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

getManifold(pp::PartialPrior{<:PackedManifoldKernelDensity}) = pp.Z.manifold

getSample(cf::CalcFactor{<:PartialPrior}) = samplePoint(cf.factor.Z)


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



