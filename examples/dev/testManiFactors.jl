using DistributedFactorGraphs
using IncrementalInference
using StaticArrays
using Manifolds


##



fg = initfg()

v0 = addVariable!(fg, :x0, Pose2)
v1 = addVariable!(fg, :x1, Point2)
v2 = addVariable!(fg, :x2, Pose2)

fg.solverParams.graphinit = false

addFactor!(fg, [:x0], Prior(Normal()))
addFactor!(fg, [:x0,:x1], Prior(Normal()))
addFactor!(fg, [:x1,:x2], Prior(Normal()))
addFactor!(fg, [:x0,:x2], Prior(Normal()))

## ======================================================================================
## 
## ======================================================================================
@defVariable Pose2 SpecialEuclidean(2) IIF.default_identity(SpecialEuclidean(2))
getManifold(Pose2)
getPointType(Pose2)
getPointIdentity(Pose2)

fg = initfg()

v0 = addVariable!(fg, :x0, Pose2)
v1 = addVariable!(fg, :x1, Pose2)


mp = ManifoldPrior(SpecialEuclidean(2),ProductRepr(SA[0., 0], SA[1.0 0; 0 1]), MvNormal([1.0, 1.0, 0.0]))
p = addFactor!(fg, [:x0], mp;  graphinit=true)

mf = ManifoldFactor(SpecialEuclidean(2), MvNormal([0.1, 0.2, 0.01]))
addFactor!(fg, [Symbol("x$i"),Symbol("x$(i+1)")], CircularCircular(Normal(1.0, 0.1))), 0:3)

f1 = addFactor!(fg, [:x0,:x1])


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

# @enter addFactor!(fg, [:x0], mp;  graphinit=true)


mf = ManifoldFactor(TranslationGroup(2), MvNormal([1.0, 2.0], [0.1,0.1]))
f = addFactor!(fg, [:x0, :x1], mf)


solveGraph!(fg)

## ======================================================================================
## 
## ======================================================================================

Base.convert(::Type{<:Tuple}, M::SpecialOrthogonal{2}) = (:Circular,)

@defVariable SO2 SpecialOrthogonal(2) [0.0]
getManifold(SO2)
getPointType(SO2)
getPointIdentity(SO2)

fg = initfg()

v0 = addVariable!(fg, :x0, SO2)
v1 = addVariable!(fg, :x1, SO2)


mp = ManifoldPrior(SpecialOrthogonal(2), SA[1.0 0.0; 0 1], MvNormal([1.0]))
p = addFactor!(fg, [:x0], mp;  graphinit=true)

doautoinit!(fg, :x0)

# @enter addFactor!(fg, [:x0], mp;  graphinit=true)


mf = ManifoldFactor(SpecialEuclidean(2), MvNormal([0.1, 0.2, 0.01]))
addFactor!(fg, [Symbol("x$i"),Symbol("x$(i+1)")], CircularCircular(Normal(1.0, 0.1))), 0:3)

f1 = addFactor!(fg, [:x0,:x1])


## ======================================================================================
## Sphere(n)
## ======================================================================================

Base.convert(::Type{<:Tuple}, M::Sphere{2, ℝ}) = (:Euclid, :Euclid)
Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{Sphere{2, ℝ}})  = (:Euclid, :Euclid)

@defVariable Sphere2 Sphere(2) [1.0, 0.0, 0.0]
M = getManifold(Sphere2)
pT = getPointType(Sphere2)
pϵ = getPointIdentity(Sphere2)

is_point(getManifold(Sphere2), getPointIdentity(Sphere2))

fg = initfg()

v0 = addVariable!(fg, :x0, Sphere2)
v1 = addVariable!(fg, :x1, Sphere2)


mp = ManifoldPrior(Sphere(2), SA[1., 0, 0], MvNormal([0.01, 0.01]))
p = addFactor!(fg, [:x0], mp)

doautoinit!(fg, :x0)
# @enter doautoinit!(fg, :x0)

# @enter addFactor!(fg, [:x0], mp;  graphinit=true)


mf = ManifoldFactor(Sphere(2), MvNormal([1.0, 2.0], [0.1,0.1]))
f = addFactor!(fg, [:x0, :x1], mf)


solveGraph!(fg)