
# using Revise
using Test
using LinearAlgebra
using Manifolds, Manopt
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

##