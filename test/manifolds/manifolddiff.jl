
# using Revise
using Test
using LinearAlgebra
using IncrementalInference
using ManifoldsBase
using Manifolds, Manopt
import Optim
using FiniteDifferences, ManifoldDiff
import Rotations as _Rot

##

# finitediff setup
r_backend = ManifoldDiff.TangentDiffBackend(
    ManifoldDiff.FiniteDifferencesBackend()
)

##
@testset "ManifoldDiff, Basic test" begin
##

# problem setup
n = 100
σ = π / 8
M = Manifolds.Sphere(2)
p = 1 / sqrt(2) * [1.0, 0.0, 1.0]
data = [exp(M, p,  σ * rand(M; vector_at=p)) for i in 1:n];

# objective function
f(M, p) = sum(1 / (2 * n) * distance.(Ref(M), Ref(p), data) .^ 2)
# f_(p) = f(M,p)

# non-manual: intrinsic finite differences gradient
function grad_f_FD(M,p)
  f_(p_) = f(M,p_)
  ManifoldDiff.gradient(M, f_, p, r_backend)
end
# manual gradient
# grad_f(M, p) = sum(1 / n * grad_distance.(Ref(M), data, Ref(p)));


# and solve
@time m1 = gradient_descent(M, f, grad_f_FD, data[1])

@info "Basic Manopt test" string(m1')
@test isapprox(p, m1; atol=0.15)

##
end

##

"""
    ManifoldWrapper{TM<:AbstractManifold} <: Optim.Manifold
    
Adapts Manifolds.jl manifolds for use in Optim.jl
"""
struct ManifoldWrapper{TM<:AbstractManifold} <: Optim.Manifold
    M::TM
end

function Optim.retract!(M::ManifoldWrapper, x)
    ManifoldsBase.embed_project!(M.M, x, x)
    return x
end

function Optim.project_tangent!(M::ManifoldWrapper, g, x)
    ManifoldsBase.embed_project!(M.M, g, x, g)
    return g
end

##
@testset "Optim.jl ManifoldWrapper example from mateuszbaran (copied to catch issues on future changes)" begin
##
# Example modified from: https://gist.github.com/mateuszbaran/0354c0edfb9cdf25e084a2b915816a09

# example usage of Manifolds.jl manifolds in Optim.jl
M = Manifolds.Sphere(2)
x0 = [1.0, 0.0, 0.0]
q = [0.0, 1.0, 0.0]

f(p) = 0.5 * distance(M, p, q)^2

# manual gradient
function g!(X, p)
    log!(M, X, p, q)
    X .*= -1
    println(p, X)
end

##

sol = Optim.optimize(f, g!, x0, Optim.ConjugateGradient(; manifold=ManifoldWrapper(M)))
@test isapprox([0,1,0.], sol.minimizer; atol=1e-6)


## finitediff gradient (non-manual)

function g_FD!(X,p)
  X .= ManifoldDiff.gradient(M, f, p, r_backend)
  X
end

#
x0 = [1.0, 0.0, 0.0]

sol = Optim.optimize(f, g_FD!, x0, Optim.ConjugateGradient(; manifold=ManifoldWrapper(M)))
@test isapprox([0,1,0.], sol.minimizer; atol=1e-6)

##

# x0 = [1.0, 0.0, 0.0]
# # internal ForwardDfif doesnt work out the box on Manifolds
# sol = Optim.optimize(f, x0, Optim.ConjugateGradient(; manifold=ManifoldWrapper(M)); autodiff=:forward )
# @test isapprox([0,1,0.], sol.minimizer; atol=1e-8)

##
end


@testset "Modified Manifolds.jl ManifoldWrapper <: Optim.Manifold for SpecialEuclidean(2; vectors=HybridTangentRepresentation())" begin
##

M = Manifolds.SpecialEuclidean(2; vectors=HybridTangentRepresentation())
e0 = ArrayPartition([0,0.], [1 0; 0 1.])

x0 = deepcopy(e0)
Cq = 9*ones(3)
while 1.5 < abs(Cq[3])
  @show Cq .= randn(3)
  # Cq[3] = 1.5 # breaks ConjugateGradient
end
q  = exp(M,e0,hat(M,e0,Cq))

f(p) = distance(M, p, q)^2

## finitediff gradient (non-manual)
function g_FD!(X,p)
  X .= ManifoldDiff.gradient(M, f, p, r_backend)
  X
end

## sanity check gradients

X = hat(M, e0, zeros(3))
g_FD!(X, q)
# gradient at the optimal point should be zero
@show X_ = [X.x[1][:]; X.x[2][:]]
@test isapprox(0, sum(abs.(X_)); atol=1e-8 )

# gradient not the optimal point should be non-zero
g_FD!(X, e0)
@show X_ = [X.x[1][:]; X.x[2][:]]
@test 0.01 < sum(abs.(X_))

## do optimization
x0 = deepcopy(e0)
sol = Optim.optimize(f, g_FD!, x0, Optim.ConjugateGradient(; manifold=ManifoldWrapper(M)))
Cq .= randn(3)
# Cq[
@show sol.minimizer
@test isapprox( f(sol.minimizer), 0; atol=1e-3 )
@test isapprox( 0, sum(abs.(log(M, e0, compose(M, inv(M,q), sol.minimizer)))); atol=1e-3)

##
end


@testset "Modified ManifoldsWrapper for Optim.Manifolds, SpecialEuclidean(3)" begin
##


M = Manifolds.SpecialEuclidean(3; vectors=HybridTangentRepresentation())
e0 = ArrayPartition([0,0,0.], Matrix(_Rot.RotXYZ(0,0,0.)))

x0 = deepcopy(e0)
Cq = 0.25*randn(6)
q  = exp(M,e0,hat(M,e0,Cq))

f(p) = distance(M, p, q)^2

## finitediff gradient (non-manual)
function g_FD!(X,p)
  X .= ManifoldDiff.gradient(M, f, p, r_backend)
  X
end

## sanity check gradients

X = hat(M, e0, zeros(6))
g_FD!(X, q)

@show X_ = [X.x[1][:]; X.x[2][:]]
# gradient at the optimal point should be zero
@test isapprox(0, sum(abs.(X_)); atol=1e-6 )

# gradient not the optimal point should be non-zero
g_FD!(X, e0)
@show X_ = [X.x[1][:]; X.x[2][:]]
@test 0.01 < sum(abs.(X_))

## do optimization
x0 = deepcopy(e0)
sol = Optim.optimize(f, g_FD!, x0, Optim.ConjugateGradient(; manifold=ManifoldWrapper(M)))
# Cq .= 0.5*randn(6)
# Cq[
@show sol.minimizer
@test isapprox( f(sol.minimizer), 0; atol=1e-3 )
@test isapprox( 0, sum(abs.(log(M, e0, compose(M, inv(M,q), sol.minimizer)))); atol=1e-3)


##
end


@testset "Optim.Manifolds, SpecialEuclidean(3), using IIF.optimizeManifold_FD" begin
##

M = Manifolds.SpecialEuclidean(3; vectors=HybridTangentRepresentation())
e0 = ArrayPartition([0,0,0.], Matrix(_Rot.RotXYZ(0,0,0.)))

x0 = deepcopy(e0)
Cq = 0.5*randn(6)
q  = exp(M,e0,hat(M,e0,Cq))

f(p) = distance(M, p, q)^2

sol = IncrementalInference.optimizeManifold_FD(M,f,x0)

@show sol.minimizer
@test isapprox( f(sol.minimizer), 0; atol=5e-3 )
@test isapprox( 0, sum(abs.(log(M, e0, compose(M, inv(M,q), sol.minimizer)))); atol=1e-3)


##
end