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

function getSample( cf::CalcFactor{<:MyFactor})
  #
  @warn "getSample(cf::CalcFactor{<:MyFactor},::Int) does not get hypo sub-selected FMD data"
  @show DFG.getLabel.(cf.metadata.fullvariables)
  # @assert DFG.getLabel.(fmd_[1].fullvariables) |> length < 3 "this factor is only between two variables"
  return rand(cf.factor.Z, 1)
end

function (cf::CalcFactor{<:MyFactor})(z, X1, X2)
  @assert DFG.getLabel.(cf.metadata.fullvariables) |> length < 3 "this factor is only between two variables"
  # just a linear difference to complete the test
  return X2 - (X1 + z)
end

##

@testset "test FactorMetadata is properly populated" begin

##


fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1_a, ContinuousScalar)
addVariable!(fg, :x1_b, ContinuousScalar)

f0 = addFactor!(fg, [:x0], Prior(Normal()))

# create the object and add it to the graph
mf = MyFactor( Normal(10,1) ) 

##

# this sampling might error
f1 = addFactor!(fg, [:x0;:x1_a;:x1_b], mf, multihypo=[1;1/2;1/2])

##

@test !isMultihypo(f0)
@test isMultihypo(f1)

##

meas = sampleFactor(fg, :x0x1_ax1_bf1, 10)

##

solveTree!(fg);

##


end

#