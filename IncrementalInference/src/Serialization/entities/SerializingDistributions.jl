
# TODO, add `<:` for concrete dispatch when using StringThemSamplableBeliefs
StringThemSamplableBeliefs = Union{
  <:Uniform,
  <:Normal,
  <:MvNormal,
  <:ZeroMeanDiagNormal,
  <:Categorical,
  <:DiscreteNonParametric,
  <:Rayleigh,
  <:BallTreeDensity,
  <:ManifoldKernelDensity,
  <:AliasingScalarSampler,
  <:HeatmapGridDensity,
  <:LevelSetGridNormal,
}

## TODO, TBD
# Base.@kwdef struct PackedDiscreteNonParametric <: PackedBelief
#   _type::String        = "IncrementalInference.PackedDiscreteNonParametric"
# end
