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
mutable struct PackedPartialPriorPassThrough <: AbstractPackedFactor
  Z::String # PackedHeatmapGridDensity
  partial::Vector{Int}
end


function convert( ::Union{Type{<:AbstractPackedFactor}, Type{<:PackedPartialPriorPassThrough}},
                  obj::PartialPriorPassThrough )
  #

  # TODO, PackedSamplableBelief
  str = convert(String, obj.Z)
  # po = convert(PackedSamplableBelief, obj.Z)
  # str = JSON2.write(po)
  PackedPartialPriorPassThrough(str, Int[obj.partial...])
end


function convert( ::Union{Type{<:AbstractFactor}, Type{<:PartialPriorPassThrough}},
                  obj::PackedPartialPriorPassThrough )
  #

  # get as bland obj to extract type
  dens = convert(SamplableBelief, obj.Z)
  # _up = JSON2.read(obj.Z)
  # _typ = DFG.getTypeFromSerializationModule(_up._type)
  # # now redo with type
  # pdens = JSON2.read(obj.Z, _typ)
  # dens = convert(SamplableBelief, pdens)
  PartialPriorPassThrough(dens, tuple(obj.partial...))
end


#