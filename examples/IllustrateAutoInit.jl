# example to illustate how autoinit functions work, and used as a development script towards a standard unit test.

using Distributions
using IncrementalInference

import IncrementalInference: getSample

mutable struct Prior{T} <: IncrementalInference.FunctorSingleton where T <: Distribution
  z::T
  Prior{T}() where {T <: Distribution} = new{T}()
  Prior(t::T) where {T <: Distribution} = new{T}(t)
  Prior{T}(t::T) where {T <: Distribution} = new{T}(t)
end
getSample(s::Prior, N::Int=1) = (reshape(rand(s.z,N),1,N), )


mutable struct LinearOffset{T} <: IncrementalInference.FunctorPairwise where T <: Distribution
  z::T
  LinearOffset{T}() where {T <: Distribution} = new{T}()
  LinearOffset(t::T) where {T <: Distribution} = new{T}(t)
  LinearOffset{T}(t::T) where {T <: Distribution} = new{T}(t)
end
getSample(s::LinearOffset, N::Int=1) = (reshape(rand(s.z,N),1,N), )

function (s::LinearOffset)(res::Array{Float64},
      idx::Int,
      meas::Tuple{Array{Float64,2}},
      X1::Array{Float64,2},
      X2::Array{Float64,2}  )
  #
  res[1] = (X2[1,idx] - X1[1,idx]) - meas[1][1,idx]
  nothing
end



fg = emptyFactorGraph()

v0 = addNode!(fg, :x0, ContinuousScalar, labels=["POSE"])

# this is unary (prior) factor and should not immediately trigger autoinit.
f1  = addFactor!(fg, [:x0], Prior(Normal())) # autoinit true throws RowVector error at setValKDE in doautoinit

# yet another uninitialized variable node
v0 = addNode!(fg, :x1, ContinuousScalar, labels=["POSE"])


lo = LinearOffset(Normal(10.0,1))

getVal(fg, :x0)

# This should call the autoinitialization procedure for :x0 and skip autoinit for :x1 given just one Pairwise factor
f1  = addFactor!(fg, [:x0, :x1], lo)


getVal(fg, :x0)
getVal(fg, :x1)

isInitialized(fg, :x0)
isInitialized(fg, :x1)


ensureAllInitialized!(fg)


using RoMEPlotting

plotKDE(fg, :x0)



isInitialized(fg, :x0)
isInitialized(fg, :x1)


plotKDE(fg, :x1)



0

## should complete and add to RoMEPlotting
# import KernelDensityEstimatePlotting: plot
# import Gadfly: plot
# import Graphs: plot
# import RoMEPlotting: plot
# function plot(fgl::FactorGraph, sym::Symbol; api::DataLayerAPI=IncrementalInference.dlapi)
#   PX = getKDE(getVert(fgl, sym, api=api))
#   plot(PX)
# end




# using Graphs
# Graphs.plot(fg.g)





#
