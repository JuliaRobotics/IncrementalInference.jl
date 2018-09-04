# development testing for multithreaded convolution

using Distributions
using IncrementalInference
using Base: Test

import IncrementalInference: getSample

struct Prior{T} <: IncrementalInference.FunctorSingleton where T <: Distribution
  z::T
end
getSample(s::Prior, N::Int=1) = (reshape(rand(s.z,N),1,:), )
struct LinearOffset{T} <: IncrementalInference.FunctorPairwise where T <: Distribution
  z::T
end
getSample(s::LinearOffset, N::Int=1) = (reshape(rand(s.z,N),1,:), )
function (s::LinearOffset)(res::Array{Float64},
                           userdata::FactorMetadata,
                           idx::Int,
                           meas::Tuple{Array{Float64, 2}},
                           X1::Array{Float64,2},
                           X2::Array{Float64,2}  )
  #
  res[1] = meas[1][idx] - (X2[1,idx] - X1[1,idx])
  nothing
end
struct MultiModalOffset <: IncrementalInference.FunctorPairwise
  z::Vector{Distribution}
  c::Categorical
end
getSample(s::MultiModalOffset, N::Int=1) = (reshape.(rand.(s.z, N),1,:)..., rand(s.c, N))
function (s::MultiModalOffset)(res::Array{Float64},
                               userdata::FactorMetadata,
                               idx::Int,
                               meas::Tuple,
                               X1::Array{Float64,2},
                               X2::Array{Float64,2}  )
  #
  res[1] = meas[meas[end][idx]][idx] - (X2[1,idx] - X1[1,idx])
  nothing
end


@testset "Basic ContinuousScalar example to ensure multithreaded convolutions work..." begin

@show Threads.nthreads()

N = 100

# Start with an empty factor graph
fg = emptyFactorGraph()

# add the first node
addNode!(fg, :x0, ContinuousScalar, N=N)

# this is unary (prior) factor and does not immediately trigger autoinit of :x0.
addFactor!(fg, [:x0], Prior(Normal(0,1)))


addNode!(fg, :x1, ContinuousScalar, N=N)
# P(Z | :x1 - :x0 ) where Z ~ Normal(10,1)
addFactor!(fg, [:x0, :x1], LinearOffset(Normal(10.0,1)))


pts = approxConv(fg, :x0x1f1, :x1, N=N)

@test 0.95*N <= sum(abs.(pts - 10.0) .< 5.0)


end


#
