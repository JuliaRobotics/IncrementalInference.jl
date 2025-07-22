
# test Manifolds

using Manifolds
using Test

##

@testset "Basic Manifolds consistency check" begin
##

w = [-0.0;-0.78;-0.18]

M = SpecialEuclideanGroup(3; variant=:right)
Mr = SpecialOrthogonalGroup(3)
pPq = ArrayPartition(zeros(3), exp(Mr, hat(LieAlgebra(Mr), w)))
rPc_ = exp(M, hat(LieAlgebra(M), [zeros(3); w], ArrayPartition{Float64}))
rPc = ArrayPartition(rPc_.x[1], rPc_.x[2])

@test isapprox(pPq.x[1], rPc.x[1])
@test isapprox(pPq.x[2], rPc.x[2])

##
end


##