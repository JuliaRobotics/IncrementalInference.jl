# prior factor that passes a density belief straight through to inference without resampling

export PartialPriorPassThrough, PackedPartialPriorPassThrough


struct PartialPriorPassThrough{B <: HeatmapDensityRegular, T <:Tuple} <: AbstractPrior
  Z::B
  partial::T
end

getManifold(pppt::PartialPriorPassThrough{<:HeatmapDensityRegular{T,H,<:ManifoldKernelDensity}}) where {T,H} = (pppt.Z.densityFnc.manifold)

# this step is skipped during main inference process
getSample(cf::CalcFactor{<:PartialPriorPassThrough}) = sampleTangent(getManifold(cf.factor), cf.factor.Z)


## ====================================================================================================
## Serialize PartialPriorPassThrough
## ====================================================================================================


mutable struct PackedPartialPriorPassThrough <: AbstractPackedFactor
  Z::PackedHeatmapDensityRegular
  partial::Vector{Int}
end


function convert( ::Union{Type{<:AbstractPackedFactor}, Type{<:PackedPartialPriorPassThrough}},
                  obj::PartialPriorPassThrough )
  #
  
  st = convert(PackedHeatmapDensityRegular, obj.Z)
  
  PackedPartialPriorPassThrough(st, Int[obj.partial...])
end


function convert( ::Union{Type{<:AbstractFactor}, Type{<:PartialPriorPassThrough}},
                  obj::PackedPartialPriorPassThrough )
  #
  
  st = convert(HeatmapDensityRegular, obj.Z)
  
  PartialPriorPassThrough(st, tuple(obj.partial...))
end


#