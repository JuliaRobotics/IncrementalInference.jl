# test new sampling API
# TODO, this test might be obsolete given later changes and requirements added to getSample interface.
#  TAC is to figure out if this functionality is still required and to remove of explain the code and related test like this file.

using Test
using IncrementalInference
using DistributedFactorGraphs

import IncrementalInference: getSample, getManifold

##

struct SpecialPrior{T <: SamplableBelief} <: AbstractPriorObservation
  z::T
end
getManifold(::SpecialPrior) = TranslationGroup(1)
getSample(s::CalcFactor{<:SpecialPrior}) = rand(s.factor.z,1)

struct SpecialLinearOffset{T <: SamplableBelief} <: RelativeObservation
  z::T
end
getManifold(::SpecialLinearOffset) = TranslationGroup(1)

function getSample(s::CalcFactor{<:SpecialLinearOffset})
  return rand(s.factor.z,1)
end


function (s::CalcFactor{<:SpecialLinearOffset})(meas,
                                                x1,
                                                x2  )
  #
  meas .- (x2 .- x1)
end

##

@testset "test specialSampler functionality..." begin

##

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0], SpecialPrior(Normal()))

addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0;:x1], SpecialLinearOffset(Normal(10,1)))

tree = solveTree!(fg)

@test calcMeanMaxSuggested(fg, :x0, :default).suggested[1] |> abs < 1.0
@test calcMeanMaxSuggested(fg, :x1, :default).suggested[1] - 10 |> abs < 3.0




## special test for IIF #568


## Singleton (Prior)

fcm = map(x->x[1], IIF._getCCW(fg, :x0f1).measurement |> deepcopy)
pts = approxConv(fg, :x0f1, :x1)
fcm2 = map(x->x[1], IIF._getCCW(fg, :x0f1).measurement)
fcm3 = map(x->x[1], IIF._getCCW(fg, :x0f1).measurement |> deepcopy)

@test 0.1 < norm(fcm - fcm2)
@test norm(fcm2 - fcm3) < 1e-5

## Pairwise

# forward direction
fcm = map(x->x[1], IIF._getCCW(fg, :x0x1f1).measurement |> deepcopy)
pts = approxConv(fg, :x0x1f1, :x1)
fcm2 = map(x->x[1], IIF._getCCW(fg, :x0x1f1).measurement)

@test 0.1 < norm(fcm - fcm2)


# reverse direction
fcm, = map(x->x[1], IIF._getCCW(fg, :x0x1f1).measurement |> deepcopy)
pts = approxConv(fg, :x0x1f1, :x0)
fcm2, = map(x->x[1], IIF._getCCW(fg, :x0x1f1).measurement)

@test 0.04 < norm(fcm - fcm2)

##

end


#
