# test null hypothesis cases

using Test
using IncrementalInference
using TensorCast
# going to introduce two new constraint types

##

@testset "Post 237 without nullhypo" begin

##

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)

addFactor!(fg, [:x0;], Prior(Normal()))
addFactor!(fg, [:x0;:x1], LinearRelative(Normal(10,1)))

solveTree!(fg)

pts_ = getBelief(fg, :x1) |> getPoints
@cast pts[i,j] := pts_[j][i]

@test 80 < sum(2 .< pts .< 18)

##

end

@testset "Post 237 nullhypo on prior" begin

##

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(10,1)), nullhypo=0.5)

solveTree!(fg)

pts_ = getBelief(fg, :x0) |> getPoints
@cast pts[i,j] := pts_[j][i]

@test 10 < sum(-15 .< pts .< 4) < 60
@test 30 < sum(4 .< pts .< 16) < 85
# @test 10 < sum(16 .< pts .< 60) < 60

##

end

@testset "Post 237 with nullhypo test" begin

##

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)

addFactor!(fg, [:x0;], Prior(Normal()))
addFactor!(fg, [:x0;:x1], LinearRelative(Normal(10,1)), nullhypo=0.5)

pts_ = approxConv(fg, :x0x1f1, :x1)
@cast pts[i,j] := pts_[j][i]

@test 20 < sum(pts .< 5)
@test 20 < sum(5 .< pts .< 15)

solveTree!(fg)

pts2_ = getBelief(fg, :x1) |> getPoints
@cast pts2[i,j] := pts2_[j][i]

@test 30 < sum(8 .< pts2 .< 12) <= 80
@test 80 < sum(-10 .< pts2 .< 30)


##

end




#
