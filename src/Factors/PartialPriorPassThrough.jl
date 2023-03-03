# prior factor that passes a density belief straight through to inference without resampling

export PartialPriorPassThrough, PackedPartialPriorPassThrough

struct PartialPriorPassThrough{
  B <: Union{<:HeatmapGridDensity, <:LevelSetGridNormal},
  T <: Tuple,
} <: AbstractPrior
  Z::B
  partial::T
end

getManifold(pppt::PartialPriorPassThrough) = getManifold(pppt.Z)

# this step is skipped during main inference process
function getSample(cf::CalcFactor{<:PartialPriorPassThrough})
  # TODO should be samplePoint for priors?
  return sampleTangent(cf.manifold, cf.factor.Z)
end

## ====================================================================================================
## Serialize PartialPriorPassThrough
## ====================================================================================================

"""
    $TYPEDEF

Required internal density to store its type
"""
Base.@kwdef mutable struct PackedPartialPriorPassThrough <: AbstractPackedFactor
  Z::PackedSamplableBelief # PackedHeatmapGridDensity
  partial::Vector{Int}
end

# StructTypes.StructType(::Type{PackedPartialPriorPassThrough}) = StructTypes.UnorderedStruct()
# StructTypes.idproperty(::Type{PackedPartialPriorPassThrough}) = :id
# StructTypes.omitempties(::Type{PackedPartialPriorPassThrough}) = (:id,)


function convert(
  ::Union{Type{<:AbstractPackedFactor}, Type{<:PackedPartialPriorPassThrough}},
  obj::PartialPriorPassThrough,
)
  #

  po = convert(PackedSamplableBelief, obj.Z)
  return PackedPartialPriorPassThrough(po, Int[obj.partial...])
end

function convert(
  ::Union{Type{<:AbstractFactor}, Type{<:PartialPriorPassThrough}},
  obj::PackedPartialPriorPassThrough,
)
  #

  dens = convert(SamplableBelief, obj.Z)
  return PartialPriorPassThrough(dens, tuple(obj.partial...))
end

#