# test if fmd in getSample and factor eval is working right



using IncrementalInference
using Test

import IncrementalInference: getSample

##

@testset "test FactorMetadata is properly populated" begin

##

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1_a, ContinuousScalar)
addVariable!(fg, :x1_b, ContinuousScalar)

addFactor!(fg, [:x0], Prior(Normal()))


struct MyFactor{T <: SamplableBelief} <: IIF.AbstractRelativeRoots
  Z::T
  # specialSampler approach will be deprecated
  specialSampler::Function
end

function getSample(mf::MyFactor, N::Int=1, fmd_...)
  @assert DFG.getLabel.(fmd_[1].fullvariables) |> length < 3 "this factor is only between two variables"
  return (reshape(rand(mf.Z, N),1,N),)
end


function (mf_::MyFactor)(res,fmd,i,meas, X1,X2)
  @assert DFG.getLabel.(fmd.fullvariables) |> length < 3 "this factor is only between two variables"

  # just a linear difference to complete the test
  res .= X2[:,i] - (X1[:,i] + meas[1][:,1])
end



##

# create the object and add it to the graph
mf = MyFactor(Normal(10,1),getSample) 

# this sampling might error
@test_broken addFactor!(fg, [:x0;:x1_a;:x1_b], mf, multihypo=[1;1/2;1/2])


@error "add factor evaluation test to testMultihypoFMD"

initManual!(fg, :x1_b, -10*ones(1,100))

# approxConv(fg, :x0f1, :x1_a)


end

#