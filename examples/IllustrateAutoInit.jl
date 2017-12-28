# example to illustate how autoinit functions work, and used as a development script towards a standard unit test.

using IncrementalInference
using Distributions

import IncrementalInference: getSample

mutable struct Prior{T} <: IncrementalInference.FunctorSingleton where T <: Distribution
  z::T
  Prior{T}() where {T <: Distribution} = new{T}()
  Prior(t::T) where {T <: Distribution} = new{T}(t)
  Prior{T}(t::T) where {T <: Distribution} = new{T}(t)
end
getSample(s::Prior, N::Int=1) = (rand(s.z,N)', )


mutable struct LinearOffset{T} <: IncrementalInference.FunctorPairwise where T <: Distribution
  z::T
  LinearOffset{T}() where {T <: Distribution} = new{T}()
  LinearOffset(t::T) where {T <: Distribution} = new{T}(t)
  LinearOffset{T}(t::T) where {T <: Distribution} = new{T}(t)
end
getSample(s::LinearOffset, N::Int=1) = (rand(s.z,N)', )

function (s::LinearOffset)(res::Array{Float64},
      idx::Int,
      meas::Tuple{Array{Float64,2}},
      X::Array{Float64,2},
      XX::Array{Float64,2}  )
  #
  res[1] = (XX[1,idx] - X[1,idx]) - meas[1][1,idx]
  nothing
end



fg = emptyFactorGraph()

# N=100

# doors = reshape(Float64[-100.0;0.0;100.0;300.0],1,4)
# pd = kde!(doors,[3.0])
# pd = resample(pd,N);
# bws = getBW(pd)[:,1]
# doors2 = getPoints(pd);

v0 = addNode!(fg, :x0, ContinuousScalar, labels=["POSE"])

ls(fg)
fg.g.vertices[1].attributes["data"]

# this is unary (prior) factor and should not immediately trigger autoinit.
f1  = addFactor!(fg, [:x0], Prior(Normal()))

# yet another uninitialized variable node
v0 = addNode!(fg, :x1, ContinuousScalar, labels=["POSE"])


lo = LinearOffset(Normal(10.0,1))

# This should call the autoinitialization procedure for :x0 and skip autoinit for :x1 given just one Pairwise factor
f1  = addFactor!(fg, [:x0, :x1], lo)







#
