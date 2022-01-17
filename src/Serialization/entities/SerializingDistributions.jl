
# TODO, add `<:` for concrete dispatch when using StringThemSamplableBeliefs
StringThemSamplableBeliefs = Union{<:Uniform, <:Normal, <:MvNormal, <:ZeroMeanDiagNormal, <:Categorical, <:DiscreteNonParametric, <:BallTreeDensity, <:ManifoldKernelDensity, <:AliasingScalarSampler, <:HeatmapGridDensity, <:LevelSetGridNormal}

## TODO, TBD
# Base.@kwdef struct PackedDiscreteNonParametric <: PackedSamplableBelief
#   _type::String        = "IncrementalInference.PackedDiscreteNonParametric"
# end

Base.@kwdef struct PackedCategorical <: PackedSamplableBelief
  _type::String        = "IncrementalInference.PackedCategorical"
  p::Vector{Float64}   = [1.0;]
end


Base.@kwdef mutable struct PackedUniform <: PackedSamplableBelief
  _type::String  = "IncrementalInference.PackedUniform"
  a::Float64     = 0.0
  b::Float64     = 1.0
  PackedSamplableTypeJSON::String = "IncrementalInference.PackedUniform"
end


Base.@kwdef struct PackedNormal <: PackedSamplableBelief
  _type::String  = "IncrementalInference.PackedNormal"
  mu::Float64    = 0.0
  sigma::Float64 = 1.0
end


Base.@kwdef struct PackedZeroMeanDiagNormal <: PackedSamplableBelief
  _type::String        = "IncrementalInference.PackedZeroMeanDiagNormal"
  diag::Vector{Float64}= ones(1)
end


Base.@kwdef struct PackedZeroMeanFullNormal <: PackedSamplableBelief
  _type::String        = "IncrementalInference.PackedZeroMeanFullNormal"
  cov::Vector{Float64} = ones(1)
end


Base.@kwdef mutable struct PackedDiagNormal <: PackedSamplableBelief
  _type::String        = "IncrementalInference.PackedDiagNormal"
  mu::Vector{Float64}  = zeros(1)
  diag::Vector{Float64}= ones(1)
end


Base.@kwdef struct PackedFullNormal <: PackedSamplableBelief
  _type::String        = "IncrementalInference.PackedFullNormal"
  mu::Vector{Float64}  = zeros(1)
  cov::Vector{Float64} = ones(1)
end


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

#