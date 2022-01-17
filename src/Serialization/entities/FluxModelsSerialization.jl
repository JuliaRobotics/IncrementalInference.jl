
export PackedFluxModelsDistribution


mutable struct PackedFluxModelsDistribution <: IIF.PackedSamplableBelief
  # standardized _type field
  _type::String
  # shape of the input data
  inputDim::Vector{Int}
  # shape of the output data
  outputDim::Vector{Int}
  # actual Flux models (Base64 encoded binary)
  mimeTypeModel::String
  models::Vector{String}
  # the data used for prediction, must be <: AbstractArray
  mimeTypeData::String
  data::String
  # shuffle model predictions relative to particle index at each sampling
  shuffle::Bool
  # false for default serialization with model info, set true for separate storage of models 
  serializeHollow::Bool
    # TODO remove requirement and standardize sampler API
    # specialSampler::Symbol
  # field name usage to direct the IIF serialization towards JSON method
  PackedSamplableTypeJSON::String
end
