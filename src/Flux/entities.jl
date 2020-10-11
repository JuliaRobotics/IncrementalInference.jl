# entities immediately available as private members in IIF, but require Flux for actual use


struct FluxModelsDistribution{ID,OD,P,D<:AbstractArray}
  # shape of the input data
  inputDim::NTuple{ID,Int}
  # shape of the output data
  outputDim::NTuple{OD,Int}
  # actual Flux models
  models::Vector{P}
  # the data used for prediction, must be <: AbstractArray
  data::D
  # shuffle model predictions relative to particle index at each sampling
  shuffle::Base.RefValue{Bool}
  # false for default serialization with model info, set true for separate storage of models 
  serializeHollow::Base.RefValue{Bool}
  # TODO remove requirement and standardize sampler API
  specialSampler::Function
end
