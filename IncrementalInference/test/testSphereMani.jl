using DistributedFactorGraphs
using IncrementalInference
using Manifolds
using StaticArrays
using Test

import Manifolds: identity_element

##

@testset "Test Sphere(2) prior and relative" begin
##

# NOTE Sphere{2} is not a lie group so the identity element does not exits.
# this is for testing only and will be removed once upgraded to support any Riemannian Manifold.
DFG.getPointIdentity(::typeof(Sphere(2))) = SVector(1.0, 0.0, 0.0)
#FIXME REMOVE! this is type piracy and not a good idea, for testing only!!!
Manifolds.identity_element(::typeof(Sphere(2))) = SVector(1.0, 0.0, 0.0)
Manifolds.identity_element(::typeof(Sphere(2)), p::AbstractVector) = SVector(1.0, 0.0, 0.0) # Float64[1,0,0]

Base.convert(::Type{<:Tuple}, M::typeof(Sphere(2))) = (:Euclid, :Euclid)
Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{typeof(Sphere(2))})  = (:Euclid, :Euclid)

@defVariable Sphere2 Sphere(2) SVector(1.0, 0.0, 0.0)
M = getManifold(Sphere2)
@test M == Sphere(2)
pT = getPointType(Sphere2)
@test pT == SVector{3,Float64}
pϵ = getPointIdentity(Sphere2)
@test pϵ == [1.0, 0.0, 0.0]

@test is_point(getManifold(Sphere2), getPointIdentity(Sphere2))

fg = initfg()

v0 = addVariable!(fg, :x0, Sphere2)

mp = ManifoldPrior(Sphere(2), SA[1., 0, 0], MvNormal(Diagonal(map(abs2, [0.01, 0.01]))), DefaultOrthonormalBasis(), ExponentialRetraction())
p = addFactor!(fg, [:x0], mp)

doautoinit!(fg, :x0)

vnd = getVariableState(fg, :x0, :default)
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

vnd = getVariableState(fg, :x1, :default)
mn_ = mean(M, vnd.val)
@info "isapprox" q mn_
@test all(isapprox.(mn_, q, atol=0.05))
@test all(is_point.(Ref(M), vnd.val))

##
end