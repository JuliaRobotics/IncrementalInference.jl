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
@show size(pts)
@test norm(Base.mean(pts,2)[1]-[0.0]) < 0.3

pts = evalFactor2(fg, f2, v1.index)
@show size(pts)
@test norm(Base.mean(pts,2)[1]-[2.0]) < 0.3


tree = wipeBuildNewTree!(fg)

inferOverTreeR!(fg,tree, N=N)
pts = getVal(fg, :x1)
@test norm(Base.mean(pts,2)[1]-[1.0]) < 0.5
@test norm(Base.mean(pts,2)[2]-[0.0]) < 0.3

# plotKDE(getVertKDE(fg, :x1),levels=3)











#
