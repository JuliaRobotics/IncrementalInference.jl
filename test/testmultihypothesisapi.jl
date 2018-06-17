using Base: Test
using IncrementalInference, Distributions
# going to introduce two new constraint types
import Base: convert
import IncrementalInference: getSample


mutable struct DevelopPrior <: FunctorSingleton
  x::Distribution
end
getSample(dpl::DevelopPrior, N::Int=1) = (reshape(rand(dpl.x, N),1,N), )

mutable struct DevelopLikelihood <: FunctorPairwise
  x::Distribution
end
getSample(dpl::DevelopLikelihood, N::Int=1) = (reshape(rand(dpl.x, N),1,N), )
function (vv::DevelopLikelihood)(res::Array{Float64},
            userdata,
            idx::Int,
            meas::Tuple,
            wXi::Array{Float64,2},
            wXj::Array{Float64,2}  )::Void
  #
  res[1] = meas[1][idx] - (wXj[1,idx] - wXi[1,idx])
  nothing
end


N  = 100
fg = emptyFactorGraph()

v1 = addNode!(fg, :x1, ContinuousScalar, N=N)

pr = DevelopPrior(Normal(10.0,1.0))
f1 = addFactor!(fg,[:x1],pr)


ensureAllInitialized!(fg)

# Juno.breakpoint("/home/dehann/.julia/v0.5/IncrementalInference/src/ApproxConv.jl",121)

pts = evalFactor2(fg, f1, v1.index, N=N)

@test sum(abs.(pts - 1.0) .< 5) < 30
@test sum(abs.(pts - 10.0) .< 5) > 30



v2 = addNode!(fg, :x2, ContinuousScalar, N=N)
pp = DevelopLikelihood(Normal(100.0,1.0))
f2 = addFactor!(fg, [:x1;:x2], pp)


ensureAllInitialized!(fg)


@test abs(Base.mean(getVal(fg, :x2))-110.0) < 10.0









v3 = addNode!(fg, :x3, ContinuousScalar, N=N)
v4 = addNode!(fg, :x4, ContinuousScalar, N=N)

ppMH = DevelopLikelihood(Normal(90.0,1.0))
f3 = addFactor!(fg, [:x2;:x3;:x4], ppMH, multihypo=(0.0,0.5,0.5))


# @test getData(f3).fnc.hypoverts == [:x3, :x4]
@test sum(abs.(getData(f3).fnc.hypotheses.p[1] .- 0.0)) < 0.1
@test sum(abs.(getData(f3).fnc.hypotheses.p[2:3] .- 0.5)) < 0.1



setVal!(v2,1*ones(1,100))
setVal!(v3,2*ones(1,100))
setVal!(v4,3*ones(1,100))


@testset "Test multi-hypothesis factor convolution exploration" begin

  pts = evalFactor2(fg, f3, v2.index, N=N)

  @test 99 < sum(pts .<= -70.0)

  pts = evalFactor2(fg, f3, v3.index, N=N)

  @test 25 < sum(pts .== 3.0) < 75

  pts = evalFactor2(fg, f3, v4.index, N=N)

  @test 25 < sum(pts .== 2.0) < 75

end



println("Packing converters")


mutable struct PackedDevelopPrior <: PackedInferenceType
  x::String
  PackedDevelopPrior() = new()
  PackedDevelopPrior(x) = new(x)
end
function convert(::Type{PackedDevelopPrior}, d::DevelopPrior)
  PackedDevelopPrior(string(d.x))
end
function convert(::Type{DevelopPrior}, d::PackedDevelopPrior)
  DevelopPrior(IncrementalInference.normalfromstring(d.x))
end

mutable struct PackedDevelopLikelihood <: PackedInferenceType
  x::String
  PackedDevelopLikelihood() = new()
  PackedDevelopLikelihood(x) = new(x)
end
function convert(::Type{PackedDevelopLikelihood}, d::DevelopLikelihood)
  PackedDevelopLikelihood(string(d.x))
end
function convert(::Type{DevelopLikelihood}, d::PackedDevelopLikelihood)
  DevelopLikelihood(IncrementalInference.extractdistribution(d.x))
end


@testset "test packing and unpacking the data structure" begin

  topack = getData(f1)
  dd = convert(PackedFunctionNodeData{PackedDevelopPrior},topack)
  unpacked = convert(FunctionNodeData{GenericWrapParam{DevelopPrior}},dd)

  @test abs(unpacked.fnc.usrfnc!.x.μ - 10.0) < 1e-10
  @test abs(unpacked.fnc.usrfnc!.x.σ - 1.0) < 1e-10



  topack = getData(f3)
  dd = convert(PackedFunctionNodeData{PackedDevelopLikelihood},topack)
  unpacked = convert(FunctionNodeData{GenericWrapParam{DevelopLikelihood}},dd)

  # @test unpacked.fnc.hypoverts == Symbol[:x3; :x4]
  @test sum(abs.(unpacked.fnc.hypotheses.p[1] .- 0.0)) < 0.1
  @test sum(abs.(unpacked.fnc.hypotheses.p[2:3] .- 0.5)) < 0.1
  # str = "Symbol[:x3, :x4];[0.5, 0.5]"
  # IncrementalInference.parsemultihypostr(str)

end

# using KernelDensityEstimatePlotting
#
#
#
# plot(kde!(pts))



#
