# partial constraint development
using Base: Test
using IncrementalInference, Distributions
# going to introduce two new constraint types
import IncrementalInference: getSample


type DevelopPartial <: IncrementalInference.FunctorSingleton
  x::Distribution
  partial::Tuple
end
getSample(dpl::DevelopPartial, N::Int=1) = (rand(dpl.x, N)', )


type DevelopDim2 <: IncrementalInference.FunctorSingleton
  x::Distribution
end
getSample(dpl::DevelopDim2, N::Int=1) = (rand(dpl.x, N), )






N=100
fg = emptyFactorGraph()


v1 = addNode!(fg,:x1,randn(2,N),N=N)

pr = DevelopDim2(MvNormal([0.0;0.0],eye(2)))
f1  = addFactor!(fg,[:x1],pr)

dp = DevelopPartial(Normal(2.0, 1.0),(1,))
f2  = addFactor!(fg,[:x1],dp)


pts = evalFactor2(fg, f1, v1.index)
# @show size(pts)
@test norm(Base.mean(pts,2)[1]-[0.0]) < 0.3

pts = evalFactor2(fg, f2, v1.index)
# @show size(pts)
@test norm(Base.mean(pts,2)[1]-[2.0]) < 0.3


tree = wipeBuildNewTree!(fg)

inferOverTreeR!(fg,tree, N=N)
pts = getVal(fg, :x1)
@test norm(Base.mean(pts,2)[1]-[1.0]) < 0.5
@test norm(Base.mean(pts,2)[2]-[0.0]) < 0.3

# plotKDE(getVertKDE(fg, :x1),levels=3)



type DevelopPartialPairwise <: IncrementalInference.FunctorPairwise
  x::Distribution
  partial::Tuple
  DevelopPartialPairwise(x::Distribution) = new(x, (2,))
end
getSample(dpl::DevelopPartialPairwise, N::Int=1) = (rand(dpl.x, N)', )

function (dp::DevelopPartialPairwise)(res::Vector{Float64},
                              idx::Int,
                              meas::Tuple{Array{Float64,2}},
                              x1::Array{Float64},
                              x2::Array{Float64}  )
  #
  res[1] = meas[1][1,idx] - (x2[2,idx]-x1[2,idx])
  nothing
end



v2 = addNode!(fg,:x2,100*randn(2,N),N=N)

valx2 = getVal(fg, :x2)

dpp = DevelopPartialPairwise(Normal(10.0, 1.0))
f3  = addFactor!(fg,[:x1;:x2],dpp)


dp2 = DevelopPartial( Normal(-20.0, 1.0), (1,) )
f4  = addFactor!(fg,[:x2],dp2)


pts = evalFactor2(fg, f3, v2.index)
plotKDE(kde!(pts))
@test size(pts,1) == 2
@test norm(Base.mean(pts,2)[2]-[10.0]) < 3.0
@test norm(valx2[1,:] - pts[1,:]) < 1e-5

pts = evalFactor2(fg, f4, v2.index)
@test size(pts,1) == 1
@test norm(Base.mean(pts,2)[1]-[-20.0]) < 0.3
@test (Base.std(pts,2)[1]-1.0) < 0.3


println("test belief prediction...")
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
@test norm(Base.mean(pts,2)[1]-[1.0]) < 0.5
@test norm(Base.mean(pts,2)[2]-[0.0]) < 0.5

pts = getVal(fg, :x2)
@test norm(Base.mean(pts,2)[1]-[-20.0]) < 2.0
@test norm(Base.mean(pts,2)[2]-[10.0]) < 2.0
@test (Base.std(pts,2)[1]-1.0) < 3.0
@test (Base.std(pts,2)[2]-1.0) < 3.0

# plotKDE(getVertKDE(fg, :x2),levels=3)

# spyCliqMat(tree, :x2)

#
