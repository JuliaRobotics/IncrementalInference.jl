# test if fmd in getSample and factor eval is working right



using IncrementalInference
using Test
import IncrementalInference: getSample, getManifold

##

struct MyFactor{T <: SamplableBelief} <: IIF.AbstractManifoldMinimize
  Z::T
  # specialSampler approach will be deprecated
  # specialSampler::Function
end

getManifold(mf::MyFactor) = TranslationGroup(getDimension(mf.Z))

function getSample( cf::CalcFactor{<:MyFactor})
  #
  @warn "getSample(cf::CalcFactor{<:MyFactor},::Int) does not get hypo sub-selected FMD data: $(DFG.getLabel.(cf.fullvariables))" cf.solvefor maxlog=1
  # @assert length( DFG.getLabel.(fmd_[1].fullvariables) ) < 3 "this factor is only between two variables"
  return rand(cf.factor.Z, 1)
end

function (cf::CalcFactor{<:MyFactor})(z, X1, X2)
  @assert length(cf.fullvariables) < 3 "this factor is only between two variables. solvefor=$(cf.solvefor)"
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

## test #424

@test_throws AssertionError addFactor!(fg, [:x0;:x1_a;:x1_b], mf, multihypo=[1/2;1/2])

##


# this sampling might error
f1 = addFactor!(fg, [:x0;:x1_a;:x1_b], mf, multihypo=[1;1/2;1/2])

##

@test !isMultihypo(f0)
@test isMultihypo(f1)

##

meas = sampleFactor(fg, :x0x1_ax1_bf1, 10)

# initAll!(fg)
# pts = approxConv(fg, :x0x1_ax1_bf1, :x1_a)
# pts = approxConv(fg, :x0x1_ax1_bf1, :x1_b)
# pts = approxConv(fg, :x0x1_ax1_bf1, :x0)

##

solveTree!(fg);

##


end

#