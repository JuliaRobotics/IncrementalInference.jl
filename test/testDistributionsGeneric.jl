
using Test
using IncrementalInference

##

@testset "test generic functions on distributions" begin
##

@test getDimension(Normal()) == 1
@test getDimension(MvNormal([1;1;0.1])) == 3


##
end