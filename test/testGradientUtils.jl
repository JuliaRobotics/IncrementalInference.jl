
using IncrementalInference
using Test

##

@testset "test _evalFactorTemporary" begin

## test utility to buidl a temporary graph

fct = EuclidDistance(Normal(10,1))
T_pt_vec = [(ContinuousScalar,[0;]); (ContinuousScalar,[9.5;])]

##

dfg, _dfgfct = IIF._buildGraphByFactorAndTypes!(fct, T_pt_vec...)


## test  the evaluation of factor without

B = IIF._evalFactorTemporary!(EuclidDistance(Normal(10,1)), 2, ([10;],), T_pt_vec... );

@test_broken B isa Vector{Vector{Float64}}
@test isapprox( B[1], [10.0;], atol=1e-6)

##

end

#