
Base.@kwdef struct PackedManifoldKernelDensity <: PackedSamplableBelief
  _type::String = "IncrementalInference.PackedManifoldKernelDensity"
  varType::String
  pts::Vector{Vector{Float64}}
  bw::Vector{Float64} = Float64[]
  partial::Vector{Int} = Int[]
  infoPerCoord::Vector{Float64} = zeros(length(pts[1]))
end

Base.@kwdef struct PackedAliasingScalarSampler <: PackedSamplableBelief
  _type::String = "IncrementalInference.PackedAliasingScalarSampler"
  domain::Vector{Float64} = [0; 1.0]
  weights::Vector{Float64} = [0.5; 0.5]
end
