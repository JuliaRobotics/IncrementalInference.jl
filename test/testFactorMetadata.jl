using IncrementalInference
using Distributions
using Base: Test



@testset "test default userdata::FactorMetadata..." begin


fgt = emptyFactorGraph()

addNode!(fgt, :x1, ContinuousScalar)
addFactor!(fgt, [:x1], Prior(Normal()))

addNode!(fgt, :x2, ContinuousScalar)
addFactor!(fgt, [:x1;:x2], LinearConditional(Normal(10,1)))

fc = getData(getVert(fgt, :x1x2f1, nt=:fnc))
@test length(fc.fnc.cpt[1].factormetadata.variableuserdata) == 2
@test fc.fnc.cpt[1].factormetadata.solvefor == :null



end







#
