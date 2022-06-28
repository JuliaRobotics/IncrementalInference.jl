using Test
using IncrementalInference

@testset "test Manual Init - distribution" begin

fg = initfg()
addVariable!(fg, :x0, ContinuousScalar)

belief = Normal(1.,0.1)
initVariable!(fg, :x0, belief)
pts = getPoints(fg, :x0)
M = getManifold(fg, :x0)
@test isapprox(mean(M, pts),[1],atol=0.1)
@test isapprox(std(M, pts),0.1,atol=0.1)

# test var api
v = getVariable(fg, :x0)
initVariable!(v, belief)
pts = getVal(v)
@test isapprox(mean(M, pts),[1],atol=0.1)
@test isapprox(std(M, pts),0.1,atol=0.1)

end