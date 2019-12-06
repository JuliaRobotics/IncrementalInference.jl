using IncrementalInference
using Test

@testset "test the basics" begin

fg = initfg()

addVariable!(fg, :x1, ContinuousScalar)
addVariable!(fg, :x2, ContinuousScalar)
addFactor!(fg, [:x1;:x2], LinearConditional(Normal()), autoinit=false)
addFactor!(fg, [:x2], Prior(Normal()), autoinit=false)

@test exists(fg, :x1)

@test !exists(fg, :l13)

end
