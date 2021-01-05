# test if fmd in getSample and factor eval is working right



using IncrementalInference
using Test

import IncrementalInference: getSample

##



struct MyFactor{T <: SamplableBelief} <: IIF.AbstractRelativeRoots
  Z::T
  # specialSampler approach will be deprecated
  # specialSampler::Function
end


function getSample(cf::CalcFactor{<:MyFactor}, N::Int=1)
  @warn "getSample(cf::CalcFactor{<:MyFactor},::Int) does not get hypo sub-selected FMD data"
  @show DFG.getLabel.(cf.metadata.fullvariables)
  # @assert DFG.getLabel.(fmd_[1].fullvariables) |> length < 3 "this factor is only between two variables"
  return (reshape(rand(cf.factor.Z, N),1,N),)
end


function (cf::CalcFactor{<:MyFactor})(res, z, X1, X2)
  @assert DFG.getLabel.(cf.metadata.fullvariables) |> length < 3 "this factor is only between two variables"

  # just a linear difference to complete the test
  res .= X2 - (X1 + z)
  nothing
end

##



@testset "test FactorMetadata is properly populated" begin


##


fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1_a, ContinuousScalar)
addVariable!(fg, :x1_b, ContinuousScalar)

addFactor!(fg, [:x0], Prior(Normal()))

# create the object and add it to the graph
mf = MyFactor( Normal(10,1) ) 

##

# this sampling might error
addFactor!(fg, [:x0;:x1_a;:x1_b], mf, multihypo=[1;1/2;1/2])

##

meas = freshSamples(fg, :x0x1_ax1_bf1, 10)

##

solveTree!(fg);

##


end

#