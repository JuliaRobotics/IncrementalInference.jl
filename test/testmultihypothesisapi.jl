# using Revise

using Test

using IncrementalInference, Distributions
using DistributedFactorGraphs
using Statistics
# going to introduce two new constraint types
import Base: convert
import IncrementalInference: getSample

##

mutable struct DevelopPrior{T <: SamplableBelief} <: AbstractPrior
  x::T
end
getSample(cf::CalcFactor{<:DevelopPrior}, N::Int=1) = (reshape(rand(cf.factor.x, N),1,N), )

mutable struct DevelopLikelihood{T <: SamplableBelief} <: AbstractRelativeRoots
  x::T
end

getSample(cf::CalcFactor{<:DevelopLikelihood}, N::Int=1) = (reshape(rand(cf.factor.x, N),1,N), )
function (cf::CalcFactor{<:DevelopLikelihood})(meas, wXi, wXj)
  #
  return meas - (wXj - wXi)
end


N  = 100
fg = initfg()

##

@testset "test populate factor graph with a multi-hypothesis factor..." begin

##

v1 = addVariable!(fg, :x1, ContinuousScalar, N=N)

pr = DevelopPrior(Normal(10.0,1.0))
f1 = addFactor!(fg,[:x1],pr)


ensureAllInitialized!(fg)

pts = evalFactor(fg, f1, v1.label, N=N)

@test sum(abs.(pts .- 1.0) .< 5) < 30
@test sum(abs.(pts .- 10.0) .< 5) > 30


v2 = addVariable!(fg, :x2, ContinuousScalar, N=N)
pp = DevelopLikelihood(Normal(100.0,1.0))
f2 = addFactor!(fg, [:x1;:x2], pp)

ensureAllInitialized!(fg)

@test abs(Statistics.mean(getVal(fg, :x2))-110.0) < 10.0


v3 = addVariable!(fg, :x3, ContinuousScalar, N=N)
v4 = addVariable!(fg, :x4, ContinuousScalar, N=N)

ppMH = DevelopLikelihood(Normal(90.0,1.0))
f3 = addFactor!(fg, [:x2;:x3;:x4], ppMH, multihypo=[1.0;0.5;0.5])


# @test IIIF._getCCW(f3).hypoverts == [:x3, :x4]
@test sum(abs.(IIF._getCCW(f3).hypotheses.p[1] .- 0.0)) < 0.1  # 1.0 becomes 0.0 for computational convenience
@test sum(abs.(IIF._getCCW(f3).hypotheses.p[2:3] .- 0.5)) < 0.1


initManual!(fg, :x2, 1*ones(1,100))
initManual!(fg, :x3, 2*ones(1,100))
initManual!(fg, :x4, 3*ones(1,100))

##

end



@testset "Test multi-hypothesis factor convolution exploration" begin

##

pts = approxConv(fg, :x2x3x4f1, :x2, N=N)

@test 99 < sum(pts .<= -70.0)

pts = approxConv(fg, :x2x3x4f1, :x3, N=N)

@test 15 < sum(pts .== 3.0) < 75

pts = approxConv(fg, :x2x3x4f1, :x4, N=N)

@test 15 < sum(pts .== 2.0) < 75

##

end

##

println("Packing converters")


mutable struct PackedDevelopPrior <: PackedInferenceType
  x::String
  PackedDevelopPrior() = new()
  PackedDevelopPrior(x) = new(x)
end
function convert(::Type{PackedDevelopPrior}, d::DevelopPrior)
  PackedDevelopPrior(convert(PackedSamplableBelief, d.x))
end
function convert(::Type{DevelopPrior}, d::PackedDevelopPrior)
  DevelopPrior(convert(SamplableBelief, d.x))
end

mutable struct PackedDevelopLikelihood <: PackedInferenceType
  x::String
  PackedDevelopLikelihood() = new()
  PackedDevelopLikelihood(x) = new(x)
end
function convert(::Type{PackedDevelopLikelihood}, d::DevelopLikelihood)
  PackedDevelopLikelihood(convert(PackedSamplableBelief, d.x))
end
function convert(::Type{DevelopLikelihood}, d::PackedDevelopLikelihood)
  DevelopLikelihood(convert(SamplableBelief, d.x))
end


##

@testset "test packing and unpacking the data structure" begin

##

topack = getSolverData(getFactor(fg,:x1f1))
dd = convert(PackedFunctionNodeData{PackedDevelopPrior},topack)
unpacked = convert(FunctionNodeData{CommonConvWrapper{DevelopPrior}},dd)

@test abs(IIF._getCCW(unpacked).usrfnc!.x.μ - 10.0) < 1e-10
@test abs(IIF._getCCW(unpacked).usrfnc!.x.σ - 1.0) < 1e-10


fct = getFactor(fg, :x2x3x4f1)
@show typeof(fct)
topack = getSolverData(fct) # f3
dd = convert(PackedFunctionNodeData{PackedDevelopLikelihood},topack)
unpacked = convert(FunctionNodeData{CommonConvWrapper{DevelopLikelihood}},dd)

# @test IIF._getCCW(unpacked).hypoverts == Symbol[:x3; :x4]
@test sum(abs.(IIF._getCCW(unpacked).hypotheses.p[1] .- 0.0)) < 0.1
@test sum(abs.(IIF._getCCW(unpacked).hypotheses.p[2:3] .- 0.5)) < 0.1


##

end

##

# start a new factor graph
N = 200
fg = initfg()

##

@testset "test tri-modal factor..." begin

##

v1 = addVariable!(fg, :x1, ContinuousScalar, N=N)

pr = DevelopPrior(Normal(10.0,1.0))
f1 = addFactor!(fg,[:x1],pr)


ensureAllInitialized!(fg)

# Juno.breakpoint("/home/dehann/.julia/v0.5/IncrementalInference/src/ApproxConv.jl",121)

pts = approxConv(fg, Symbol(f1.label), :x1, N=N)


@test sum(abs.(pts .- 1.0) .< 5) < 30
@test sum(abs.(pts .- 10.0) .< 5) > 30



v2 = addVariable!(fg, :x2, ContinuousScalar, N=N)
pp = DevelopLikelihood(Normal(100.0,1.0))
f2 = addFactor!(fg, [:x1;:x2], pp)

ensureAllInitialized!(fg)

@test abs(Statistics.mean(getVal(fg, :x2))-110.0) < 10.0



v3 = addVariable!(fg, :x3, ContinuousScalar, N=N)
v4 = addVariable!(fg, :x4, ContinuousScalar, N=N)
v5 = addVariable!(fg, :x5, ContinuousScalar, N=N)


ppMH = DevelopLikelihood(Normal(90.0,1.0))
f3 = addFactor!(fg, [:x2;:x3;:x4;:x5], ppMH, multihypo=[1.0,0.333,0.333,0.334])



# @test IIF._getCCW(f3).hypoverts == [:x3, :x4]
@test sum(abs.(IIF._getCCW(f3).hypotheses.p[1] .- 0.0)) < 0.1  # 1.0 becomes 0.0 for computational convenience
@test sum(abs.(IIF._getCCW(f3).hypotheses.p[2] .- 0.333)) < 0.001
@test sum(abs.(IIF._getCCW(f3).hypotheses.p[3] .- 0.333)) < 0.001
@test sum(abs.(IIF._getCCW(f3).hypotheses.p[4] .- 0.334)) < 0.001


initManual!(fg, :x2 ,1*ones(1,100))
initManual!(fg, :x3 ,2*ones(1,100))
initManual!(fg, :x4 ,3*ones(1,100))
initManual!(fg, :x5 ,4*ones(1,100))


# solve for certain idx
pts = approxConv(fg, :x2x3x4x5f1, :x2, N=N)

@test 0.95*N < sum(pts .<= -70.0)


# solve for one of uncertain variables
pts = approxConv(fg, :x2x3x4x5f1, :x3, N=N)

@test 0.1*N < sum(80 .< pts .< 100.0) < 0.5*N
@test 0.1*N < sum(pts .== 3.0) < 0.5*N
@test 0.1*N < sum(pts .== 4.0) < 0.5*N

# 0.7 to accomodate bad-init null hypo
@test 0.5*N <= sum(70 .< pts .< 110.0) + sum(pts .== 3.0) + sum(pts .== 4.0)


# solve for one of uncertain variables
pts = approxConv(fg, :x2x3x4x5f1, :x4, N=N)

@test 0.1*N < sum(80 .< pts .< 100.0) < 0.5*N
@test 0.1*N < sum(pts .== 2.0) < 0.5*N
@test 0.1*N < sum(pts .== 4.0) < 0.5*N

@test 0.5*N <= sum(80 .< pts .< 100.0) + sum(pts .== 2.0) + sum(pts .== 4.0)


# solve for one of uncertain variables
pts = approxConv(fg, :x2x3x4x5f1, :x5, N=N)

@test 0.1*N < sum(80 .< pts .< 100.0) < 0.5*N
@test 0.1*N < sum(pts .== 2.0) < 0.5*N
@test 0.1*N < sum(pts .== 3.0) < 0.5*N

@test 0.5*N <= sum(80 .< pts .< 100.0) + sum(pts .== 2.0) + sum(pts .== 3.0)

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
