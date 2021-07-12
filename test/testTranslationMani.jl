using DistributedFactorGraphs
using IncrementalInference
using Manifolds
using StaticArrays
using Test

@testset "Test TranslationGroup(2)" begin
## ======================================================================================
## 
## ======================================================================================

Base.convert(::Type{<:Tuple}, M::TranslationGroup{Tuple{2},ℝ}) = (:Euclid, :Euclid)
Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{TranslationGroup{Tuple{2},ℝ}})  = (:Euclid, :Euclid)

@defVariable Point2 TranslationGroup(2) [0.0, 0.0]
getManifold(Point2)
getPointType(Point2)
getPointIdentity(Point2)

fg = initfg()

v0 = addVariable!(fg, :x0, Point2)
v1 = addVariable!(fg, :x1, Point2)


mp = ManifoldPrior(TranslationGroup(2), SA[10., 20], MvNormal([1.0,1.0]))
p = addFactor!(fg, [:x0], mp;  graphinit=true)

doautoinit!(fg, :x0)

mf = ManifoldFactor(TranslationGroup(2), MvNormal([1.0, 2.0], [0.1,0.1]))
f = addFactor!(fg, [:x0, :x1], mf)

solveGraph!(fg)

end