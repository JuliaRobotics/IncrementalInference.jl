
# using Revise
using Test
using LinearAlgebra
using ManifoldsBase
using Manifolds, Manopt
using Optim
using FiniteDifferences, ManifoldDiff

##

@testset "ManifoldDiff, Basic test" begin
##

# finitediff setup
r_backend = ManifoldDiff.TangentDiffBackend(
    ManifoldDiff.FiniteDifferencesBackend()
)

# problem setup
n = 100
σ = π / 8
M = Sphere(2)
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


@testset "Optim.jl ManifoldWrapper example from mateuszbaran (copied to catch issues on future changes)" begin
# Example modified from: https://gist.github.com/mateuszbaran/0354c0edfb9cdf25e084a2b915816a09
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

sol = optimize(f, g!, x0, ConjugateGradient(; manifold=ManifoldWrapper(M)))
@test isapprox([0,1,0.], sol.minimizer; atol=1e-8)


## finitediff gradient (non-manual)

r_backend = ManifoldDiff.TangentDiffBackend(
    ManifoldDiff.FiniteDifferencesBackend()
)
function g_FD!(X,p)
  X .= ManifoldDiff.gradient(M, f, p, r_backend)
  X
end

#
x0 = [1.0, 0.0, 0.0]

sol = optimize(f, g_FD!, x0, ConjugateGradient(; manifold=ManifoldWrapper(M)))
@test isapprox([0,1,0.], sol.minimizer; atol=1e-8)

##

# x0 = [1.0, 0.0, 0.0]
# # internal ForwardDfif doesnt work out the box on Manifolds
# sol = optimize(f, x0, ConjugateGradient(; manifold=ManifoldWrapper(M)); autodiff=:forward )
# @test isapprox([0,1,0.], sol.minimizer; atol=1e-8)

##
end

##