
Base.@kwdef struct PackedCategorical <: PackedBelief
  _type::String = "IncrementalInferenceTypes.PackedCategorical"
  p::Vector{Float64} = [1.0;]
end

Base.@kwdef mutable struct PackedUniform <: PackedBelief
  _type::String = "IncrementalInferenceTypes.PackedUniform"
  a::Float64 = 0.0
  b::Float64 = 1.0
  PackedSamplableTypeJSON::String = "IncrementalInferenceTypes.PackedUniform"
end

Base.@kwdef struct PackedNormal <: PackedBelief
  _type::String = "IncrementalInferenceTypes.PackedNormal"
  mu::Float64 = 0.0
  sigma::Float64 = 1.0
end

Base.@kwdef struct PackedZeroMeanDiagNormal <: PackedBelief
  _type::String = "IncrementalInferenceTypes.PackedZeroMeanDiagNormal"
  diag::Vector{Float64} = ones(1)
end

Base.@kwdef struct PackedZeroMeanFullNormal <: PackedBelief
  _type::String = "IncrementalInferenceTypes.PackedZeroMeanFullNormal"
  cov::Vector{Float64} = ones(1)
end

Base.@kwdef mutable struct PackedDiagNormal <: PackedBelief
  _type::String = "IncrementalInferenceTypes.PackedDiagNormal"
  mu::Vector{Float64} = zeros(1)
  diag::Vector{Float64} = ones(1)
end

Base.@kwdef struct PackedFullNormal <: PackedBelief
  _type::String = "IncrementalInferenceTypes.PackedFullNormal"
  mu::Vector{Float64} = zeros(1)
  cov::Vector{Float64} = ones(1)
end

Base.@kwdef struct PackedRayleigh <: PackedBelief
  _type::String = "IncrementalInferenceTypes.PackedRayleigh"
  sigma::Float64 = 1.0
end

#
