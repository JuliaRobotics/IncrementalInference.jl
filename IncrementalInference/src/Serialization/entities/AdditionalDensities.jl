
Base.@kwdef struct PackedManifoldKernelDensity <: PackedBelief
  _type::String = "IncrementalInference.PackedManifoldKernelDensity"
  varType::String
  pts::Vector{Vector{Float64}}
  bw::Vector{Float64} = Float64[]
  partial::Vector{Int} = Int[]
  infoPerCoord::Vector{Float64} = zeros(length(pts[1]))
end

Base.@kwdef struct PackedAliasingScalarSampler <: PackedBelief
  _type::String = "IncrementalInference.PackedAliasingScalarSampler"
  domain::Vector{Float64} = [0; 1.0]
  weights::Vector{Float64} = [0.5; 0.5]
end


function DFG.pack(dtr::AliasingScalarSampler)
  return PackedAliasingScalarSampler(; domain = dtr.domain, weights = dtr.weights.values)
end

function DFG.unpack(dtr::PackedAliasingScalarSampler)
  return AliasingScalarSampler(dtr.domain, dtr.weights ./ sum(dtr.weights))
end
