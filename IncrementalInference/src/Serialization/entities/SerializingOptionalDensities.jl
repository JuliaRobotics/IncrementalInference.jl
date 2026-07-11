


Base.@kwdef mutable struct PackedHeatmapGridDensity <: PackedBelief
  _type::String = "IncrementalInference.PackedHeatmapGridDensity"
  data::Vector{Vector{Float64}}
  domain::Tuple{Vector{Float64}, Vector{Float64}}
  hint_callback::String
  bw_factor::Float64
  N::Int
  # _densityFnc::String = "" # only use if storing parched belief data entry label/id
end


Base.@kwdef mutable struct PackedLevelSetGridNormal <: PackedBelief
  _type::String = "IncrementalInference.PackedLevelSetGridNormal"
  level::Float64
  sigma::Float64
  sigma_scale::Float64
  # make sure the JSON nested packing works with the serialization overlords
  heatmap::PackedHeatmapGridDensity
end


Base.@kwdef mutable struct PackedFluxModelsDistribution <: PackedBelief
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
  # TODO, only use ._type.  Legacy, field name usage to direct the IIF serialization towards JSON method
  PackedSamplableTypeJSON::String
end
