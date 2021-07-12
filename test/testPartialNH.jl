# test nullhypo with n-dim partial

using Test
using IncrementalInference

##

@testset "test nullhypo with n-dimensional partial" begin

##

fg = initfg()
addVariable!(fg, :x0, ContinuousEuclid{3})

addFactor!(fg, [:x0;], PartialPrior(MvNormal(zeros(2), ones(2)),(2,3)), nullhypo=0.2)

addVariable!(fg, :x1, ContinuousEuclid{3})
addFactor!(fg, [:x1;], PartialPrior(Normal(10,1),(1,)))
addFactor!(fg, [:x0; :x1], LinearRelative(MvNormal([10;0;0.0], ones(3))), nullhypo=0.2)

##

solveTree!(fg);

##

@warn "Suppressing testPartialNH.jl during transition to Manifolds.jl"
@test isapprox( getPPE(fg, :x0).suggested, [0;0;0], atol=1)
@test isapprox( getPPE(fg, :x1).suggested, [10;0;0], atol=1)

##

end

#