using DistributedFactorGraphs
using IncrementalInference
using Manifolds
using StaticArrays
using Test

##

@testset "Test Sphere(2) prior and relative" begin
##

# NOTE us getPointIdentity instead
#FIXME REMOVE! this is type piracy and not a good idea, for testing only!!!
# Manifolds.identity_element(::Sphere{2, ℝ}, p::Vector{Float64}) = Float64[1,0,0]

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

mp = ManifoldPrior(Sphere(2), SA[1., 0, 0], MvNormal(Diagonal(map(abs2, [0.01, 0.01]))), DefaultOrthonormalBasis(), ExponentialRetraction())
p = addFactor!(fg, [:x0], mp)

doautoinit!(fg, :x0)

vnd = getVariableSolverData(fg, :x0)
@test all(isapprox.(mean(M, vnd.val), [1,0,0], atol=0.1))
@test all(is_point.(Ref(M), vnd.val))

v1 = addVariable!(fg, :x1, Sphere2)
mf = ManifoldFactor(Sphere(2), MvNormal([0.1, 0.2], [0.05,0.05]))
f = addFactor!(fg, [:x0, :x1], mf)

##

smtasks = Task[]
solveTree!(fg; smtasks)

#
p = SA[1.,0,0]
X = get_vector(M, p, SA[0.1,0.2], DefaultOrthonormalBasis())
q = exp(M, p, X)

vnd = getVariableSolverData(fg, :x1)
@test all(isapprox.(mean(M, vnd.val), q, atol=0.01))
@test all(is_point.(Ref(M), vnd.val))

##
end