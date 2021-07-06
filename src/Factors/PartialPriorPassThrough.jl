# prior factor that passes a density belief straight through to inference without resampling

export PartialPriorPassThrough


struct PartialPriorPassThrough{B <: HeatmapDensityRegular, T <:Tuple} <: AbstractPrior
  Z::B
  partial::T
end

# this step is skipped during main inference process
getSample(cf::CalcFactor{<:PartialPriorPassThrough}, N::Int=1) = (reshape(rand(cf.factor.Z.densityFnc,N),:,N), )