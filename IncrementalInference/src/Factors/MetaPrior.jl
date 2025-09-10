# Meta prior brings additional information not necessarily numerical prior

Base.@kwdef struct MetaPrior{T} <: AbstractPriorObservation
  data::T
  partial::Vector{Int} = Int[]
end
MetaPrior(data) = MetaPrior(;data)

getManifold(::MetaPrior) = TranslationGroup(0)
getMeasurementParametric(::MetaPrior) = MvNormal(zeros(0))

getSample(cf::CalcFactor{<:MetaPrior}) = SVector()