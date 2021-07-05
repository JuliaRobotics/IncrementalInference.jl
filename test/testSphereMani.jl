using DistributedFactorGraphs
using IncrementalInference
using Manifolds
using StaticArrays
using Test

@testset "Test Sphere(2) prior" begin

Base.convert(::Type{<:Tuple}, M::Sphere{2, ℝ}) = (:Euclid, :Euclid)
Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{Sphere{2, ℝ}})  = (:Euclid, :Euclid)

@defVariable Sphere2 Sphere(2) [1.0, 0.0, 0.0]
M = getManifold(Sphere2)
@test M == Sphere(2)
pT = getPointType(Sphere2)
@test pT == Vector{Float64}
pϵ = getPointIdentity(Sphere2)
@test pϵ == [1.0, 0.0, 0.0]

@test is_point(getManifold(Sphere2), getPointIdentity(Sphere2))

fg = initfg()

v0 = addVariable!(fg, :x0, Sphere2)

mp = ManifoldPrior(Sphere(2), SA[1., 0, 0], MvNormal([0.01, 0.01]))
p = addFactor!(fg, [:x0], mp)

doautoinit!(fg, :x0)

vnd = getVariableSolverData(fg, :x0)
@test all(isapprox.(mean(vnd.val), [1,0,0], atol=0.1))
@test all(is_point.(Ref(M), vnd.val))

v1 = addVariable!(fg, :x1, Sphere2)
mf = ManifoldFactor(Sphere(2), MvNormal([0.1, 0.2], [0.05,0.05]))
f = addFactor!(fg, [:x0, :x1], mf)

##
# Debugging Sphere error
# DimensionMismatch("tried to assign 200 elements to 300 destinations")
smtasks = Task[]
solveTree!(fg; smtasks, verbose=true, recordcliqs=ls(fg))
hists = fetchCliqHistoryAll!(smtasks);

end