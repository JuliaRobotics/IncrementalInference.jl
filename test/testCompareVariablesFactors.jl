using Revise

using IncrementalInference
using Test



@testset "testing compare functions for variables and factors..." begin

fg = initfg()

v1 = addVariable!(fg, :x0, ContinuousScalar)
f1 = addFactor!(fg, [:x0;], Prior(Normal()))

@test compareVariable(v1,v1)
@test compareFactor(f1,f1)

v2 = addVariable!(fg, :x1, ContinuousScalar)
f2 = addFactor!(fg, [:x0;:x1], LinearConditional(Normal(2.0, 0.1)))

fg2 = deepcopy(fg)


@test !compareVariable(v1,v2)

@test !compareFactor(f1,f2)

@test compareAllVariables(fg, fg)

@test compareAllVariables(fg, fg2)

@test compareSimilarVariables(fg, fg)
@test compareSimilarVariables(fg, fg2)

@test compareSimilarFactors(fg, fg)
@test compareSimilarFactors(fg, fg2)

@test compareFactorGraphs(fg, fg)
@test compareFactorGraphs(fg, fg2)

end







#
