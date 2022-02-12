
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



#