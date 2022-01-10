
# TODO, add `<:` for concrete dispatch when using StringThemSamplableBeliefs
StringThemSamplableBeliefs = Union{Normal, MvNormal, ZeroMeanDiagNormal, Categorical, DiscreteNonParametric, BallTreeDensity, ManifoldKernelDensity, AliasingScalarSampler}


Base.@kwdef struct PackedCategorical <: PackedSamplableBelief
  _type::String        = "IncrementalInference.PackedCategorical"
end


Base.@kwdef struct PackedDiscreteNonParametric <: PackedSamplableBelief
  _type::String        = "IncrementalInference.PackedDiscreteNonParametric"
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
  cov::Matrix{Float64} = ones(1,1)
end


Base.@kwdef mutable struct PackedDiagNormal <: PackedSamplableBelief
  _type::String        = "IncrementalInference.PackedMvNormal"
  mu::Vector{Float64}  = [0.0;]
  diag::Vector{Float64}= ones(1)
end


Base.@kwdef struct PackedFullNormal <: PackedSamplableBelief
  _type::String        = "IncrementalInference.PackedFullNormal"
  mu::Vector{Float64}  = [0.0;]
  cov::Matrix{Float64} = ones(1,1)
end


Base.@kwdef struct PackedAliasingScalarSampler <: PackedSamplableBelief
  _type::String        = "IncrementalInference.PackedAliasingScalarSampler"
end




#