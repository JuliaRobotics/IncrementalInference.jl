


Base.@kwdef struct PackedManifoldKernelDensity <: PackedSamplableBelief
  _type::String            = "IncrementalInference.PackedManifoldKernelDensity"
  varType::String
  pts::Vector{Vector{Float64}}
  bw::Vector{Float64}      = ones(length(pts[1]))
  partial::Vector{Int}     = collect(1:length(pts[1])) 
  infoPerCoord::Vector{Float64} = zeros(length(pts[1]))
end


Base.@kwdef struct PackedAliasingScalarSampler <: PackedSamplableBelief
  _type::String           = "IncrementalInference.PackedAliasingScalarSampler"
  domain::Vector{Float64} = [0;1.0;]
  weights::Vector{Float64}= [0.5;0.5;]
end


Base.@kwdef mutable struct PackedHeatmapGridDensity <: PackedSamplableBelief
  _type::String = "IncrementalInference.PackedHeatmapGridDensity"
  data::Vector{Vector{Float64}}
  domain::Tuple{Vector{Float64}, Vector{Float64}}
  hint_callback::String
  bw_factor::Float64
  N::Int
  # _densityFnc::String = "" # only use if storing parched belief data entry label/id
end


Base.@kwdef mutable struct PackedLevelSetGridNormal <: PackedSamplableBelief
  _type::String = "IncrementalInference.PackedLevelSetGridNormal"
  level::Float64
  sigma::Float64
  sigma_scale::Float64
  # make sure the JSON nested packing works with the serialization overlords
  heatmap::PackedHeatmapGridDensity
end