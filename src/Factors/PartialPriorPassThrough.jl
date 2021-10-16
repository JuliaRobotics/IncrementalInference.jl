# prior factor that passes a density belief straight through to inference without resampling

export PartialPriorPassThrough, PackedPartialPriorPassThrough


struct PartialPriorPassThrough{B <: Union{<:HeatmapGridDensity,<:LevelSetGridNormal}, T <:Tuple} <: AbstractPrior
  Z::B
  partial::T
end

getManifold(pppt::PartialPriorPassThrough) = getManifold(pppt.Z)
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

  str = JSON2.write(obj.Z)
  PackedPartialPriorPassThrough(str, Int[obj.partial...])
end


function convert( ::Union{Type{<:AbstractFactor}, Type{<:PartialPriorPassThrough}},
                  obj::PackedPartialPriorPassThrough )
  #
  
  _typ = DFG.getTypeFromSerializationModule(obj.Z._type)
  dens = JSON2.read(obj.Z, _typ)
  
  PartialPriorPassThrough(dens, tuple(obj.partial...))
end


#