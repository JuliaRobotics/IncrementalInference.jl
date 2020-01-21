# partial constraint development
using Test
using IncrementalInference
# going to introduce two new constraint mutable structs
import IncrementalInference: getSample
using Statistics


mutable struct DevelopPartial <: IncrementalInference.FunctorSingleton
  x::Distribution
  partial::Tuple # should rather be static types or templates for performance Tuple{Int, Int} etc.
end
getSample(dpl::DevelopPartial, N::Int=1) = (rand(dpl.x, N)', )


mutable struct DevelopDim2 <: IncrementalInference.FunctorSingleton
  x::Distribution
end
getSample(dpl::DevelopDim2, N::Int=1) = (rand(dpl.x, N), )





N=100 # 50
fg = initfg()


v1 = addVariable!(fg,:x1,ContinuousMultivariate(2,manifolds=(:Euclid, :Euclid)),N=N)

pr = DevelopDim2(MvNormal([0.0;0.0], 0.01*Matrix{Float64}(LinearAlgebra.I,2,2)))
f1  = addFactor!(fg,[:x1],pr)

dp = DevelopPartial(Normal(2.0, 1.0),(1,))
f2  = addFactor!(fg,[:x1],dp)


@testset "test evaluation of full constraint prior" begin
  pts = evalFactor2(fg, f1, v1.label, N=N)
  @test size(pts,1) == 2
  @test size(pts,2) == N
  @test norm(Statistics.mean(pts,dims=2)[1] .- [0.0]) < 0.3
end

memcheck = getVal(v1)

@testset "test evaluation of partial constraint prior" begin
  X1pts = getVal(v1)
  pts = evalFactor2(fg, f2, v1.label, N=N)

  @test size(pts, 1) == 2
  @test size(pts, 2) == N
  @test norm(Statistics.mean(pts,dims=2)[1] .- [2.0]) < 0.75
  # ensure the correct response from
  @test norm(X1pts[1,:] - pts[1,:]) > 2.0
  @test norm(X1pts[2,:] - pts[2,:]) < 1e-10
  @test norm(X1pts - memcheck) < 1e-10
end

@testset "test solving of factor graph" begin
  tree = wipeBuildNewTree!(fg)

  inferOverTree!(fg,tree, N=N)
  pts = getVal(fg, :x1)
  @test norm(Statistics.mean(pts,dims=2)[1] .- [0.0]) < 0.25
  @test norm(Statistics.mean(pts,dims=2)[2] .- [0.0]) < 0.25
end
# plotKDE(getVertKDE(fg, :x1),levels=3)



mutable struct DevelopPartialPairwise <: IncrementalInference.FunctorPairwise
  x::Distribution
  partial::Tuple
  DevelopPartialPairwise(x::Distribution) = new(x, (2,))
end
getSample(dpl::DevelopPartialPairwise, N::Int=1) = (rand(dpl.x, N)', )

function (dp::DevelopPartialPairwise)(res::Vector{Float64},
                              userdata::FactorMetadata,
                              idx::Int,
                              meas::Tuple, #{RowVector{Float64,Array{Float64,1}}}, #Tuple{Array{Float64,2}},
                              x1::Array{Float64},
                              x2::Array{Float64}  )
  #
  res[1] = meas[1][1,idx] - (x2[2,idx]-x1[2,idx])
  nothing
end



v2 = addVariable!(fg,:x2,ContinuousMultivariate(2),N=N)


dpp = DevelopPartialPairwise(Normal(10.0, 1.0))
f3  = addFactor!(fg,[:x1;:x2],dpp)


dp2 = DevelopPartial( Normal(-20.0, 1.0), (1,) )
f4  = addFactor!(fg,[:x2;], dp2, autoinit=false)
doautoinit!(fg, :x2)

# drawGraph(fg, show=true)


@testset "test evaluation of multiple simultaneous partial constraints" begin
global fg

ensureAllInitialized!(fg)
valx2 = getVal(fg, :x2)
pts = approxConv(fg, :x1x2f1, :x2, N=N) # evalFactor2(fg, f3, v2.index, N=N)
@test size(pts,1) == 2
@test norm(Statistics.mean(pts,dims=2)[2] .- [10.0]) < 3.0
@test norm(valx2[1,:] - pts[1,:]) < 1e-5

pts = approxConv(fg, :x2f1, :x2, N=N) # evalFactor2(fg, f4, v2.index, N=N)
@test size(pts,1) == 2
@test norm(Statistics.mean(pts,dims=2)[1] .- [-20.0]) < 0.75
@test (Statistics.std(pts,dims=2)[1] .- 1.0) < 0.4

end

# keep previous values to ensure funciton evaluation is modifying correct data fields

@warn "restore findRelatedFromPotential as testset!"
# @testset "test findRelatedFromPotential..." begin
# global v2, fg, f3, f4, N

thefac = getFactor(fg, :x1x2f1)

X2lpts = getVal(getVariable(fg, :x2))
keepaside, = findRelatedFromPotential(fg, thefac, :x2, N)
@test Ndim(keepaside) == 2
lpts = KernelDensityEstimate.getPoints(keepaside)
@test size(lpts,2) == N

@show X2lpts[2,95:100]
@show lpts[2,95:100]
@show getPoints(keepaside)

# DevelopPartialPairwise must only modify the second dimension of proposal distribution on X2
@test norm(X2lpts[1,:] - lpts[1,:]) < 1e-10
# @test norm(X2lpts[2,:] - lpts[2,:]) > 1e-10 # 10*N*0.5 # why so big?
memcheck = getVal(v2)
@test norm(X2lpts - memcheck) < 1e-10


X2lpts = getVal(v2)
p4, = findRelatedFromPotential(fg, f4, v2.label, N)
@test Ndim(p4) == 2
lpts = KernelDensityEstimate.getPoints(keepaside)
@test size(lpts,2) == N

# DevelopPartialPairwise must only modify the second dimension of proposal distribution on X2
@test norm(X2lpts[1,:] - lpts[1,:]) < 1e-10
@test norm(X2lpts[2,:] - lpts[2,:]) > 1e-10 # 10*N*0.5 # why so big?
memcheck = getVal(v2)
@test norm(X2lpts - memcheck) < 1e-10

# end


@testset "test belief prediction with partials..." begin

global v2, fg

# partial prior
X2pts = getVal(v2)
val, = predictbelief(fg, v2, [f4], N=N)
@test norm(X2pts[2,:] - val[2,:]) < 1e-10
@test 0.0 < norm(X2pts[1,:] - val[1,:])
@test norm(Statistics.mean(val[1,:]) .+ 20.0) < 0.75


# partial pairwise
X2pts = getVal(v2)
val, = predictbelief(fg, v2, [f3], N=N)
@test norm(X2pts[1,:] - val[1,:]) < 1e-10
@test 0.0 < norm(X2pts[2,:] - val[2,:])
@test abs(Statistics.mean(val[2,:] - getVal(v1)[2,:]) .- 10.0) < 0.75


# combination of partials
val, = predictbelief(fg, v2, [f3;f4], N=N)
# plotKDE(kde!(val),levels=3)
@test norm(Statistics.mean(val,dims=2)[1] .- [-20.0]) < 2.0
@test norm(Statistics.mean(val,dims=2)[2] .- [10.0]) < 2.0
@test (Statistics.std(val,dims=2)[1] .- 1.0) < 3.0
@test (Statistics.std(val,dims=2)[2] .- 1.0) < 3.0


tree = wipeBuildNewTree!(fg )#, drawpdf=true)
# run(`evince bt.pdf`)

inferOverTree!(fg,tree, N=N)

pts = getVal(fg, :x1)
@test norm(Statistics.mean(pts,dims=2)[1] .- [0.0]) < 0.5
@test norm(Statistics.mean(pts,dims=2)[2] .- [0.0]) < 0.5

pts = getVal(fg, :x2)
@test norm(Statistics.mean(pts,dims=2)[1] .- [-20.0]) < 2.0
@test norm(Statistics.mean(pts,dims=2)[2] .- [10.0]) < 2.0
@test (Statistics.std(pts,dims=2)[1]-1.0) < 3.0
@test (Statistics.std(pts,dims=2)[2]-1.0) < 3.0

end

# plotKDE(getVertKDE(fg, :x2),levels=3)

# spyCliqMat(tree, :x2)

#
