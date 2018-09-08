# default node types in IIF

SamplableBelief = Union{Distributions.Distribution, KernelDensityEstimate.BallTreeDensity, AliasingScalarSampler}


struct ContinuousScalar <: InferenceVariable
  dims::Int
  labels::Vector{String}
  ContinuousScalar() = new(1, String[])
end
struct ContinuousMultivariate <:InferenceVariable
  dims::Int
  labels::Vector{String}
  ContinuousMultivariate() = new()
  ContinuousMultivariate(x) = new(x, String[])
end


struct Prior{T} <: IncrementalInference.FunctorSingleton where T <: SamplableBelief
  z::T
end
getSample(s::Prior, N::Int=1) = (reshape(rand(s.z,N),:,N), )

struct LinearConditional{T} <: IncrementalInference.FunctorPairwise where T <: SamplableBelief
  z::T
end
getSample(s::LinearConditional, N::Int=1) = (reshape(rand(s.z,N),:,N), )
function (s::LinearConditional)(res::Array{Float64},
                                userdata::FactorMetadata,
                                idx::Int,
                                meas::Tuple{Array{Float64, 2}},
                                X1::Array{Float64,2},
                                X2::Array{Float64,2}  )
  #
  res[1] = meas[1][idx] - (X2[1,idx] - X1[1,idx])
  nothing
end
