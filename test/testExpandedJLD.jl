using IncrementalInference
using Distributions
using Base: Test



@testset "test solve from loading jld..." begin


fgt = emptyFactorGraph()

addNode!(fgt, :x1, ContinuousScalar)
addFactor!(fgt, [:x1], Prior(Normal()))

addNode!(fgt, :x2, ContinuousScalar)
addFactor!(fgt, [:x1;:x2], LinearConditional(Normal(10,1)))

savejld(fgt)
fgl, = loadjld()

fct = getData(fgt, :x1x2f1, nt=:fct).fnc.cpt[1].factormetadata
fcl = getData(fgl, :x1x2f1, nt=:fct).fnc.cpt[1].factormetadata

@test fct.solvefor == fcl.solvefor
@test length(fct.variableuserdata) == length(fcl.variableuserdata)
@test fct.variableuserdata == fcl.variableuserdata

batchSolve!(fgt)
batchSolve!(fgl)

# getData(fgt, :x1).softtype

end



@testset "test backwards compatibility of loadjld from previous files..." begin



fc = getData(getVert(fgt, :x1x2f1, nt=:fnc))
@test length(fc.fnc.cpt[1].factormetadata.variableuserdata) == 2
@test fc.fnc.cpt[1].factormetadata.solvefor == :null



end








#
