# test null hypothesis cases

using Test
using IncrementalInference
# going to introduce two new constraint types
# import IncrementalInference: getSample

@testset "Post 237 without nullhypo" begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)

addFactor!(fg, [:x0;], Prior(Normal()))
addFactor!(fg, [:x0;:x1], LinearConditional(Normal(10,1)))

solveTree!(fg)

pts_ = getBelief(fg, :x1) |> getPoints

@test 80 < sum(2 .< pts_ .< 18)

end


@testset "Post 237 nullhypo on prior" begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(10,1)), nullhypo=0.5)

solveTree!(fg)

pts_ = getBelief(fg, :x0) |> getPoints

@test 10 < sum(-8 .< pts_ .< 8) < 80
@test 30 < sum(5 .< pts_ .< 15) < 85

end

@testset "Post 237 with nullhypo test" begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)

addFactor!(fg, [:x0;], Prior(Normal()))
addFactor!(fg, [:x0;:x1], LinearConditional(Normal(10,1)), nullhypo=0.5)


pts = approxConv(fg, :x0x1f1, :x1)

@test 20 < sum(pts .< 5)
@test 20 < sum(5 .< pts .< 15)

solveTree!(fg)

pts2 = getBelief(fg, :x1) |> getPoints

@test 30 < sum(8 .< pts2 .< 12) <= 80
@test 80 < sum(-10 .< pts2 .< 30)

end




#
