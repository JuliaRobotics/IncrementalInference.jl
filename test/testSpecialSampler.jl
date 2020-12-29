# test new sampling API

using Test
using IncrementalInference
using DistributedFactorGraphs

# import IncrementalInference: getSample


struct SpecialPrior{T} <: AbstractPrior where T <: Distribution
  z::T
  specialSampler::Function
end
specialSample(s::SpecialPrior, N::Int, fmd::FactorMetadata, Xi) = (reshape(rand(s.z,N),1,:), )
SpecialPrior(z::T) where {T <: SamplableBelief}= SpecialPrior{T}(z, specialSample)

struct SpecialLinearOffset{T <: SamplableBelief} <: AbstractRelativeRoots
  z::T
  specialSampler::Function
end

specialSample(s::SpecialLinearOffset, N::Int, fmd::FactorMetadata, Xi, Xj) = (reshape(rand(s.z,N),1,:), )

SpecialLinearOffset(z::T) where {T <: SamplableBelief} = SpecialLinearOffset{T}(z, specialSample)

function (s::SpecialLinearOffset)(res::Array{Float64},
                                  userdata::FactorMetadata,
                                  idx::Int,
                                  meas::Tuple{Array{Float64, 2}},
                                  X1::Array{Float64,2},
                                  X2::Array{Float64,2}  )
  #
  res[1] = meas[1][idx] - (X2[1,idx] - X1[1,idx])
  nothing
end


@testset "test specialSampler functionality..." begin

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

fcm, = getSolverData(getFactor(fg, :x0f1)).fnc.measurement |> deepcopy
pts = approxConv(fg, :x0f1, :x1)
fcm2, = getSolverData(getFactor(fg, :x0f1)).fnc.measurement
fcm3, = getSolverData(getFactor(fg, :x0f1)).fnc.measurement |> deepcopy

@test 0.1 < norm(fcm - fcm2)
@test norm(fcm2 - fcm3) < 1e-5

## Pairwise

# forward direction
fcm, = getSolverData(getFactor(fg, :x0x1f1)).fnc.measurement |> deepcopy
pts = approxConv(fg, :x0x1f1, :x1)
fcm2, = getSolverData(getFactor(fg, :x0x1f1)).fnc.measurement

@test 0.1 < norm(fcm - fcm2)


# reverse direction
fcm, = getSolverData(getFactor(fg, :x0x1f1)).fnc.measurement |> deepcopy
pts = approxConv(fg, :x0x1f1, :x0)
fcm2, = getSolverData(getFactor(fg, :x0x1f1)).fnc.measurement

@test 0.1 < norm(fcm - fcm2)


0
end


#
