using IncrementalInference
using Test



@testset "test default userdata::FactorMetadata..." begin


fgt = initfg()

addVariable!(fgt, :x1, ContinuousScalar)
addFactor!(fgt, [:x1], Prior(Normal()))

addVariable!(fgt, :x2, ContinuousScalar)
addFactor!(fgt, [:x1;:x2], LinearConditional(Normal(10,1)))

fc = solverData(getFactor(fgt, :x1x2f1))
@test length(fc.fnc.cpt[1].factormetadata.variableuserdata) == 2
@test fc.fnc.cpt[1].factormetadata.solvefor == :null


end








#
