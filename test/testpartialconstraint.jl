# partial constraint development
using Base: Test
using IncrementalInference, Distributions
# going to introduce two new constraint mutable structs
import IncrementalInference: getSample


mutable struct DevelopPartial <: IncrementalInference.FunctorSingleton
  x::Distribution
  partial::Tuple
end
getSample(dpl::DevelopPartial, N::Int=1) = (rand(dpl.x, N)', )


mutable struct DevelopDim2 <: IncrementalInference.FunctorSingleton
  x::Distribution
end
getSample(dpl::DevelopDim2, N::Int=1) = (rand(dpl.x, N), )





N=50
fg = emptyFactorGraph()


v1 = addNode!(fg,:x1,ContinuousMultivariate(2),N=N)

pr = DevelopDim2(MvNormal([0.0;0.0],0.01*eye(2)))
f1  = addFactor!(fg,[:x1],pr)

dp = DevelopPartial(Normal(2.0, 1.0),(1,))
f2  = addFactor!(fg,[:x1],dp)


@testset "test evaluation of full constraint prior" begin
  pts = evalFactor2(fg, f1, v1.index, N=N)
  @test size(pts,1) == 2
  @test size(pts,2) == N
  @test norm(Base.mean(pts,2)[1]-[0.0]) < 0.3
end

memcheck = getVal(v1)

@testset "test evaluation of partial constraint prior" begin
  X1pts = getVal(v1)
  pts = evalFactor2(fg, f2, v1.index, N=N)

  @test size(pts, 1) == 2
  @test size(pts,2) == N
  @test norm(Base.mean(pts,2)[1]-[2.0]) < 0.75
  # ensure the correct response from
  @test norm(X1pts[1,:] - pts[1,:]) > 2.0
  @test norm(X1pts[2,:] - pts[2,:]) < 1e-10
  @test norm(X1pts - memcheck) < 1e-10
end

@testset "test solving of factor graph" begin
  tree = wipeBuildNewTree!(fg)

  inferOverTreeR!(fg,tree, N=N)
  pts = getVal(fg, :x1)
  @test norm(Base.mean(pts,2)[1]-[0.0]) < 0.25
  @test norm(Base.mean(pts,2)[2]-[0.0]) < 0.25
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



v2 = addNode!(fg,:x2,ContinuousMultivariate(2),N=N)


dpp = DevelopPartialPairwise(Normal(10.0, 1.0))
f3  = addFactor!(fg,[:x1;:x2],dpp)


dp2 = DevelopPartial( Normal(-20.0, 1.0), (1,) )
f4  = addFactor!(fg,[:x2;], dp2)


@testset "test evaluation of multiple simultaneous partial constraints" begin

ensureAllInitialized!(fg)
valx2 = getVal(fg, :x2)
pts = evalFactor2(fg, f3, v2.index, N=N)
@test size(pts,1) == 2
@test norm(Base.mean(pts,2)[2]-[10.0]) < 3.0
@test norm(valx2[1,:] - pts[1,:]) < 1e-5

pts = evalFactor2(fg, f4, v2.index, N=N)
@test size(pts,1) == 2
@test norm(Base.mean(pts,2)[1]-[-20.0]) < 0.75
@test (Base.std(pts,2)[1]-1.0) < 0.4

end

# keep previous values to ensure funciton evaluation is modifying correct data fields

@testset "test findRelatedFromPotential..." begin

X2pts = getVal(v2)
p3 = findRelatedFromPotential(fg, f3, v2.index, N)
@test Ndim(p3) == 2
pts = KernelDensityEstimate.getPoints(p3)
@test size(pts,2) == N


# DevelopPartialPairwise must only modify the second dimension of proposal distribution on X2
@test norm(X2pts[1,:] - pts[1,:]) < 1e-10
@test norm(X2pts[2,:] - pts[2,:]) > 1e-10 # 10*N*0.5 # why so big?
memcheck = getVal(v2)
@test norm(X2pts - memcheck) < 1e-10


X2pts = getVal(v2)
p4 = findRelatedFromPotential(fg, f4, v2.index, N)
@test Ndim(p4) == 2
pts = KernelDensityEstimate.getPoints(p3)
@test size(pts,2) == N

# DevelopPartialPairwise must only modify the second dimension of proposal distribution on X2
@test norm(X2pts[1,:] - pts[1,:]) < 1e-10
@test norm(X2pts[2,:] - pts[2,:]) > 1e-10 # 10*N*0.5 # why so big?
memcheck = getVal(v2)
@test norm(X2pts - memcheck) < 1e-10

end


@testset "test belief prediction with partials..." begin

# partial prior
X2pts = getVal(v2)
val = predictbelief(fg, v2, [f4], N=N)
@test norm(X2pts[2,:] - val[2,:]) < 1e-10
@test 0.0 < norm(X2pts[1,:] - val[1,:])
@test norm(Base.mean(val[1,:])+20.0) < 0.75


# partial pairwise
X2pts = getVal(v2)
val = predictbelief(fg, v2, [f3], N=N)
@test norm(X2pts[1,:] - val[1,:]) < 1e-10
@test 0.0 < norm(X2pts[2,:] - val[2,:])
@test abs(Base.mean(val[2,:] - getVal(v1)[2,:])-10.0) < 0.75


# combination of partials
val = predictbelief(fg, v2, [f3;f4], N=N)
# plotKDE(kde!(val),levels=3)
@test norm(Base.mean(val,2)[1]-[-20.0]) < 2.0
@test norm(Base.mean(val,2)[2]-[10.0]) < 2.0
@test (Base.std(val,2)[1]-1.0) < 3.0
@test (Base.std(val,2)[2]-1.0) < 3.0


tree = wipeBuildNewTree!(fg )#, drawpdf=true)
# run(`evince bt.pdf`)

inferOverTreeR!(fg,tree, N=N)

pts = getVal(fg, :x1)
@test norm(Base.mean(pts,2)[1]-[0.0]) < 0.5
@test norm(Base.mean(pts,2)[2]-[0.0]) < 0.5

pts = getVal(fg, :x2)
@test norm(Base.mean(pts,2)[1]-[-20.0]) < 2.0
@test norm(Base.mean(pts,2)[2]-[10.0]) < 2.0
@test (Base.std(pts,2)[1]-1.0) < 3.0
@test (Base.std(pts,2)[2]-1.0) < 3.0

end

# plotKDE(getVertKDE(fg, :x2),levels=3)

# spyCliqMat(tree, :x2)

#
