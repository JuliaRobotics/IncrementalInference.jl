
# test Manifolds

using Manifolds
using Test

##

@testset "Basic Manifolds consistency check" begin
##

w = [-0.0;-0.78;-0.18]

M = SpecialEuclidean(3; vectors=HybridTangentRepresentation())
Mr = M.manifold[2]
pPq = ArrayPartition(zeros(3), exp(Mr, Identity(Mr), hat(Mr, Identity(Mr), w)))
rPc_ = exp(M, Identity(M), hat(M, Identity(M), [zeros(3);w]))
rPc = ArrayPartition(rPc_.x[1], rPc_.x[2])

@test isapprox(pPq.x[1], rPc.x[1])
@test isapprox(pPq.x[2], rPc.x[2])

##
end


##