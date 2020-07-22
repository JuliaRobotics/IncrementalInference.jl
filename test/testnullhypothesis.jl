# test null hypothesis cases

using Test
using IncrementalInference
# going to introduce two new constraint types
# import IncrementalInference: getSample


@testset "Post 237 nullhypo test" begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)

addFactor!(fg, [:x0;], Prior(Normal()))
addFactor!(fg, [:x0;:x1], LinearConditional(Normal(10,1)), nullhypo=0.5)


pts = approxConv(fg, :x0x1f1, :x1)

@test 20 < sum(pts .< 5)
@test 20 < sum(5 .< pts .< 15)

end




#
