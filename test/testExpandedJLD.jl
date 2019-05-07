using IncrementalInference
using Test

# setSerializationNamespace!("Main" => Main)

@testset "test solve by saving and loading basic jld..." begin


global fgt = FactorGraph()

addVariable!(fgt, :x1, ContinuousScalar)
addFactor!(fgt, [:x1], Prior(Normal()))

addVariable!(fgt, :x2, ContinuousScalar)
addFactor!(fgt, [:x1;:x2], LinearConditional(Normal(10,1)))

savejld(fgt)
# savejld(fgt, file=joinpath(Pkg.dir("IncrementalInference"),"test","testdata","compatibility_test_fg.jld"))
global fgl, = loadjld()

global fct = getData(fgt, :x1x2f1, nt=:fct).fnc.cpt[1].factormetadata
global fcl = getData(fgl, :x1x2f1, nt=:fct).fnc.cpt[1].factormetadata

@test fct.solvefor == fcl.solvefor
@test length(fct.variableuserdata) == length(fcl.variableuserdata)

batchSolve!(fgt)
batchSolve!(fgl)

# test again after solve
@test fct.solvefor == fcl.solvefor
@test length(fct.variableuserdata) == length(fcl.variableuserdata)

# save load and repeat tests
savejld(fgt, file="tempfg.jld2")
global fgl, = loadjld(file="tempfg.jld2")

Base.rm("tempfg.jld2")

global fct = getData(fgt, :x1x2f1, nt=:fct).fnc.cpt[1].factormetadata
global fcl = getData(fgl, :x1x2f1, nt=:fct).fnc.cpt[1].factormetadata

# @test fct.solvefor == fcl.solvefor # defaults to :null during load
@test length(fct.variableuserdata) == length(fcl.variableuserdata)


end



@testset "test backwards compatibility of loadjld from previous versions of IIF (from an existing test data file)..." begin

global filename = joinpath(dirname(pathof(IncrementalInference)), "..", "test","testdata","compatibility_test_fg.jld2")
# savejld(fgt, file=filename)  # when changing to a new compat test file (human required)
global fgprev, = loadjld(file=filename )

global fc = getData(getVert(fgprev, :x1x2f1, nt=:fnc))
@test length(fc.fnc.cpt[1].factormetadata.variableuserdata) == 2
@test fc.fnc.cpt[1].factormetadata.solvefor == :null

batchSolve!(fgprev)


end








#
