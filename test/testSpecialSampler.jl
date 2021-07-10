# test new sampling API

using Test
using IncrementalInference
using DistributedFactorGraphs

import IncrementalInference: getSample

##

struct SpecialPrior{T <: SamplableBelief} <: AbstractPrior
  z::T
end
getSample(s::CalcFactor{<:SpecialPrior}, N::Int) = ([rand(s.factor.z,1) for _ in 1:N], )

struct SpecialLinearOffset{T <: SamplableBelief} <: AbstractRelativeRoots
  z::T
end

function getSample(s::CalcFactor{<:SpecialLinearOffset}, N::Int)
  fmd = s.metadata
  ret = [rand(s.factor.z,1) for _ in 1:N]
  (ret, )
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

tree, smt, hist = solveTree!(fg)

@test getPPE(fg, :x0).suggested[1] |> abs < 1.0

@test getPPE(fg, :x1).suggested[1] - 10 |> abs < 3.0




## special test for IIF #568


## Singleton (Prior)

fcm, = IIF._getCCW(fg, :x0f1).measurement |> deepcopy
pts = approxConv(fg, :x0f1, :x1)
fcm2, = IIF._getCCW(fg, :x0f1).measurement
fcm3, = IIF._getCCW(fg, :x0f1).measurement |> deepcopy

@test 0.1 < norm(fcm - fcm2)
@test norm(fcm2 - fcm3) < 1e-5

## Pairwise

# forward direction
fcm, = IIF._getCCW(fg, :x0x1f1).measurement |> deepcopy
pts = approxConv(fg, :x0x1f1, :x1)
fcm2, = IIF._getCCW(fg, :x0x1f1).measurement

@test 0.1 < norm(fcm - fcm2)


# reverse direction
fcm, = IIF._getCCW(fg, :x0x1f1).measurement |> deepcopy
pts = approxConv(fg, :x0x1f1, :x0)
fcm2, = IIF._getCCW(fg, :x0x1f1).measurement

@test 0.1 < norm(fcm - fcm2)

##

end


#
