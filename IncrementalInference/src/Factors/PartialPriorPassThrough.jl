# prior factor that passes a density belief straight through to inference without resampling

export PartialPriorPassThrough#, PackedPartialPriorPassThrough

DFG.@tags struct PartialPriorPassThrough{
  B <: Union{<:HeatmapGridDensity, <:LevelSetGridNormal},
  T <: Tuple,
} <: AbstractPriorObservation
  Z::B & (lower = DFG.Packed, choosetype = DFG.resolvePackedType)
  partial::T & (choosetype = x->NTuple{length(x), Int},)
end

DFG.getManifold(pppt::PartialPriorPassThrough) = getManifold(pppt.Z)

# this step is skipped during main inference process
function getSample(cf::CalcFactor{<:PartialPriorPassThrough})
  # TODO should be samplePoint for priors?
  return sampleTangent(cf.manifold, cf.factor.Z)
end

## ====================================================================================================
## Serialize PartialPriorPassThrough
## ====================================================================================================

# """
#     $TYPEDEF

# Required internal density to store its type
# """
# Base.@kwdef mutable struct PackedPartialPriorPassThrough <: AbstractPackedObservation
#   Z::PackedBelief # PackedHeatmapGridDensity
#   partial::Vector{Int}
# end

# function convert(
#   ::Union{Type{<:AbstractPackedObservation}, Type{<:PackedPartialPriorPassThrough}},
#   obj::PartialPriorPassThrough,
# )
#   #

#   po = convert(PackedBelief, obj.Z)
#   return PackedPartialPriorPassThrough(po, Int[obj.partial...])
# end

# function convert(
#   ::Union{Type{<:AbstractObservation}, Type{<:PartialPriorPassThrough}},
#   obj::PackedPartialPriorPassThrough,
# )
#   #

#   dens = convert(SamplableBelief, obj.Z)
#   return PartialPriorPassThrough(dens, tuple(obj.partial...))
# end

#