


struct PackedManifoldKernelDensity <: PackedSamplableBelief
  _type::String
  varType::String
  pts::Vector{Vector{Float64}}
  bw::Vector{Float64}
  partial::Vector{Int}
  infoPerCoord::Vector{Float64}
end


Base.@kwdef struct PackedAliasingScalarSampler <: PackedSamplableBelief
  _type::String           = "IncrementalInference.PackedAliasingScalarSampler"
  domain::Vector{Float64} = [0;1.0;]
  weights::Vector{Float64}= [0.5;0.5;]
end


mutable struct PackedHeatmapGridDensity <: PackedSamplableBelief
  _type::String
  data::Vector{Vector{Float64}}
  domain::Tuple{Vector{Float64}, Vector{Float64}}
  hint_callback::String
  bw_factor::Float64
  N::Int
  # densityFnc::String # TODO rather rebuild at unpack
end


mutable struct PackedLevelSetGridNormal <: PackedSamplableBelief
  _type::String
  level::Float64
  sigma::Float64
  sigma_scale::Float64
  # make sure the JSON nested packing works with the serialization overlords
  heatmap::PackedHeatmapGridDensity
end