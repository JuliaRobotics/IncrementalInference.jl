# using Revise

using Test

using IncrementalInference, Distributions
using DistributedFactorGraphs
using Statistics
using TensorCast
# going to introduce two new constraint types
import Base: convert
import IncrementalInference: getSample, getManifold

##

mutable struct DevelopPrior{T <: SamplableBelief} <: AbstractPriorObservation
  # keeping to test user case using `.x` rather than default `.Z`
  x::T
end
getManifold(dp::DevelopPrior) = TranslationGroup(getDimension(dp.x))
getSample(cf::CalcFactor{<:DevelopPrior}) = rand(cf.factor.x, 1)

mutable struct DevelopLikelihood{T <: SamplableBelief} <: RelativeObservation
  x::T
end

getManifold(dp::DevelopLikelihood) = TranslationGroup(getDimension(dp.x))
getSample(cf::CalcFactor{<:DevelopLikelihood}) = rand(cf.factor.x, 1)
(cf::CalcFactor{<:DevelopLikelihood})(meas, wXi, wXj) = meas - (wXj - wXi)


##

N  = 100
fg = initfg()

##

@testset "test populate factor graph with a multi-hypothesis factor..." begin
##

v1 = addVariable!(fg, :x1, ContinuousScalar, N=N)

pr = DevelopPrior(Normal(10.0,1.0))
f1 = addFactor!(fg,[:x1],pr)


initAll!(fg)

pts_, _ = evalFactor(fg, f1, v1.label, N=N)
@cast pts[i,j] := pts_[j][i]

@test sum(abs.(pts .- 1.0) .< 5) < 30
@test sum(abs.(pts .- 10.0) .< 5) > 30


v2 = addVariable!(fg, :x2, ContinuousScalar, N=N)
pp = DevelopLikelihood(Normal(100.0,1.0))
f2 = addFactor!(fg, [:x1;:x2], pp)

initAll!(fg)

pts_ = getVal(fg, :x2)
@cast pts[i,j] := pts_[j][i]
@test abs(Statistics.mean(pts)-110.0) < 10.0


v3 = addVariable!(fg, :x3, ContinuousScalar, N=N)
v4 = addVariable!(fg, :x4, ContinuousScalar, N=N)

ppMH = DevelopLikelihood(Normal(90.0,1.0))
f3 = addFactor!(fg, [:x2;:x3;:x4], ppMH, multihypo=[1.0;0.5;0.5])


# @test IIIF._getCCW(f3).hypoverts == [:x3, :x4]
@test sum(abs.(IIF._getCCW(f3).hyporecipe.hypotheses.p[1] .- 0.0)) < 0.1  # 1.0 becomes 0.0 for computational convenience
@test sum(abs.(IIF._getCCW(f3).hyporecipe.hypotheses.p[2:3] .- 0.5)) < 0.1


initVariable!(fg, :x2, [1*ones(1) for _ in 1:N])
initVariable!(fg, :x3, [2*ones(1) for _ in 1:N])
initVariable!(fg, :x4, [3*ones(1) for _ in 1:N])

##
end


@testset "Test multi-hypothesis factor convolution exploration" begin
##

pts_ = approxConv(fg, :x2x3x4f1, :x2, N=N)
@cast pts[i,j] := pts_[j][i]
@test 99 < sum(pts .<= -70.0)

pts_ = approxConv(fg, :x2x3x4f1, :x3, N=N)
@cast pts[i,j] := pts_[j][i]
15 < sum(70 .< pts .< 110) < 75

pts_ = approxConv(fg, :x2x3x4f1, :x4, N=N)
@cast pts[i,j] := pts_[j][i]
15 < sum(70 .< pts .< 110) < 75


##

end

##

println("Packing converters")


Base.@kwdef mutable struct PackedDevelopPrior <: AbstractPackedObservation
  x::PackedBelief
end
function convert(::Type{PackedDevelopPrior}, d::DevelopPrior)
  PackedDevelopPrior(convert(PackedBelief, d.x))
end
function convert(::Type{DevelopPrior}, d::PackedDevelopPrior)
  DevelopPrior(convert(SamplableBelief, d.x))
end

@kwdef mutable struct PackedDevelopLikelihood <: AbstractPackedObservation
  x::PackedBelief
end
function convert(::Type{PackedDevelopLikelihood}, d::DevelopLikelihood)
  PackedDevelopLikelihood(convert(PackedBelief, d.x))
end
function convert(::Type{<:DevelopLikelihood}, d::PackedDevelopLikelihood)
  DevelopLikelihood(convert(SamplableBelief, d.x))
end


##

@testset "test packing and unpacking the data structure" begin

##
packedfac = packFactor(getFactor(fg,:x1f1))
unpackedfac = unpackFactor(packedfac)

rebuildFactorCache!(fg, unpackedfac)

@test abs(IIF._getCCW(unpackedfac).usrfnc!.x.μ - 10.0) < 1e-10
@test abs(IIF._getCCW(unpackedfac).usrfnc!.x.σ - 1.0) < 1e-10

#
packedfac = packFactor(getFactor(fg,:x2x3x4f1))
unpackedfac = unpackFactor(packedfac)
rebuildFactorCache!(fg, unpackedfac)

@test sum(abs.(IIF._getCCW(unpackedfac).hyporecipe.hypotheses.p[1] .- 0.0)) < 0.1
@test sum(abs.(IIF._getCCW(unpackedfac).hyporecipe.hypotheses.p[2:3] .- 0.5)) < 0.1

# fct = getFactor(fg, :x2x3x4f1)
# topack = getState(fct) # f3
# dd = convert(PackedFunctionNodeData{PackedDevelopLikelihood},topack)
# unpacked = reconstFactorData(fg, [:x2;:x3;:x4], FunctionNodeData{CommonConvWrapper{DevelopLikelihood}},dd)
# @test sum(abs.(IIF._getCCW(unpacked).hyporecipe.hypotheses.p[1] .- 0.0)) < 0.1
# @test sum(abs.(IIF._getCCW(unpacked).hyporecipe.hypotheses.p[2:3] .- 0.5)) < 0.1

##

end

##

# start a new factor graph
N = 200
fg = initfg()
getSolverParams(fg).N = N

##

@testset "test tri-modal factor..." begin

##

v1 = addVariable!(fg, :x1, ContinuousScalar, N=N)

pr = DevelopPrior(Normal(10.0,1.0))
f1 = addFactor!(fg,[:x1],pr)


initAll!(fg)

@test length(getVal(fg, :x1)) == N 

pts_ = approxConv(fg, Symbol(f1.label), :x1, N=N)
@cast pts[i,j] := pts_[j][i]

@test sum(abs.(pts .- 1.0) .< 5) < 30
@test sum(abs.(pts .- 10.0) .< 5) > 30



v2 = addVariable!(fg, :x2, ContinuousScalar, N=N)
pp = DevelopLikelihood(Normal(100.0,1.0))
f2 = addFactor!(fg, [:x1;:x2], pp)

initAll!(fg)
pts_ = getVal(fg, :x2)
@cast pts[i,j] := pts_[j][i]
@test abs(Statistics.mean(pts)-110.0) < 10.0



v3 = addVariable!(fg, :x3, ContinuousScalar, N=N)
v4 = addVariable!(fg, :x4, ContinuousScalar, N=N)
v5 = addVariable!(fg, :x5, ContinuousScalar, N=N)


ppMH = DevelopLikelihood(Normal(90.0,1.0))
f3 = addFactor!(fg, [:x2;:x3;:x4;:x5], ppMH, multihypo=[1.0,0.333,0.333,0.334])



# @test IIF._getCCW(f3).hypoverts == [:x3, :x4]
@test sum(abs.(IIF._getCCW(f3).hyporecipe.hypotheses.p[1] .- 0.0)) < 0.1  # 1.0 becomes 0.0 for computational convenience
@test sum(abs.(IIF._getCCW(f3).hyporecipe.hypotheses.p[2] .- 0.333)) < 0.001
@test sum(abs.(IIF._getCCW(f3).hyporecipe.hypotheses.p[3] .- 0.333)) < 0.001
@test sum(abs.(IIF._getCCW(f3).hyporecipe.hypotheses.p[4] .- 0.334)) < 0.001


initVariable!(fg, :x2 ,[1*ones(1) for _ in 1:N])
initVariable!(fg, :x3 ,[2*ones(1) for _ in 1:N])
initVariable!(fg, :x4 ,[3*ones(1) for _ in 1:N])
initVariable!(fg, :x5 ,[4*ones(1) for _ in 1:N])


# solve for certain idx
pts_ = approxConv(fg, :x2x3x4x5f1, :x2, N=N)
@cast pts[i,j] := pts_[j][i]
@test 0.95*N < sum(pts .<= -70.0)


# solve for one of uncertain variables
pts_ = approxConv(fg, :x2x3x4x5f1, :x3, N=N)
@cast pts[i,j] := pts_[j][i]
@test 0.1*N < sum(80 .< pts .< 100.0) < 0.5*N
# @test 0.1*N < sum(pts .== 3.0) < 0.5*N
# @test 0.1*N < sum(pts .== 4.0) < 0.5*N

# 0.7 to accomodate bad-init null hypo
# @test 0.5*N <= sum(70 .< pts .< 110.0) + sum(pts .== 3.0) + sum(pts .== 4.0)


# solve for one of uncertain variables
pts_ = approxConv(fg, :x2x3x4x5f1, :x4, N=N)
@cast pts[i,j] := pts_[j][i]
@test 0.1*N < sum(80 .< pts .< 100.0) < 0.5*N
# @test 0.1*N < sum(pts .== 2.0) < 0.5*N
# @test 0.1*N < sum(pts .== 4.0) < 0.5*N

# @test 0.5*N <= sum(80 .< pts .< 100.0) + sum(pts .== 2.0) + sum(pts .== 4.0)


# solve for one of uncertain variables
pts_ = approxConv(fg, :x2x3x4x5f1, :x5, N=N)
@cast pts[i,j] := pts_[j][i]
@test 0.1*N < sum(80 .< pts .< 100.0) < 0.5*N
# @test 0.1*N < sum(pts .== 2.0) < 0.5*N
# @test 0.1*N < sum(pts .== 3.0) < 0.5*N

# @test 0.5*N <= sum(80 .< pts .< 100.0) + sum(pts .== 2.0) + sum(pts .== 3.0)

##
end


@testset "test multihypo api numerical tolerance, #1086" begin
##

fg = initfg()

addVariable!(fg, :x0, ContinuousEuclid{1})
addVariable!(fg, :x1a, ContinuousEuclid{1})
addVariable!(fg, :x1b, ContinuousEuclid{1})
addFactor!(fg, [:x0;:x1a;:x1b], LinearRelative(Normal()), multihypo=[1; 0.5;0.4999999999999])
addFactor!(fg, [:x0;:x1a;:x1b], LinearRelative(Normal()), multihypo=[1; 0.5;0.5000000000001])


##
end



#
