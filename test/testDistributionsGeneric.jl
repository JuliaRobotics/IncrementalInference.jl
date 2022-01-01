
using Test
using IncrementalInference

##

@testset "test generic functions on distributions" begin
##

@test getDimension(Uniform(0,1)) == 1
@test getDimension(Normal()) == 1
@test getDimension(MvNormal([1;1;0.1])) == 3

p = manikde!(ContinuousScalar, [randn(1) for i in 1:100])
@test getDimension(p) == 1

##
end