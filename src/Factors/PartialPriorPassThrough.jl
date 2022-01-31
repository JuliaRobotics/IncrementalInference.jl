# prior factor that passes a density belief straight through to inference without resampling

export PartialPriorPassThrough, PackedPartialPriorPassThrough


struct PartialPriorPassThrough{B <: Union{<:HeatmapGridDensity,<:LevelSetGridNormal}, T <:Tuple} <: AbstractPrior
  Z::B
  partial::T
end

getManifold(pppt::PartialPriorPassThrough) = getManifold(pppt.Z)

# this step is skipped during main inference process
getSample(cf::CalcFactor{<:PartialPriorPassThrough}) = sampleTangent(getManifold(cf.factor), cf.factor.Z)


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


function convert( ::Union{Type{<:AbstractPackedFactor}, Type{<:PackedPartialPriorPassThrough}},
                  obj::PartialPriorPassThrough )
  #

  po = convert(PackedSamplableBelief, obj.Z)
  PackedPartialPriorPassThrough(po, Int[obj.partial...])
end


function convert( ::Union{Type{<:AbstractFactor}, Type{<:PartialPriorPassThrough}},
                  obj::PackedPartialPriorPassThrough )
  #

  dens = convert(SamplableBelief, obj.Z)
  PartialPriorPassThrough(dens, tuple(obj.partial...))
end


#