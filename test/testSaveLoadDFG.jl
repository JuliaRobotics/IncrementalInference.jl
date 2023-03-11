
using IncrementalInference
using Test

##

@testset "saving to and loading from FileDFG" begin
##

fg = generateGraph_Kaess()
addVariable!(fg, :x4, ContinuousScalar)
addFactor!(fg, [:x2;:x3;:x4], LinearRelative(Normal()), multihypo=[1.0;0.6;0.4])

saveFolder = "/tmp/dfg_test"
saveDFG(fg, saveFolder)

retDFG = initfg()
retDFG = loadDFG!(retDFG, saveFolder)
Base.rm(saveFolder*".tar.gz")

@test symdiff(ls(fg), ls(retDFG)) == []
@test symdiff(lsf(fg), lsf(retDFG)) == []

@show getFactor(fg, :x2x3x4f1).solverData.multihypo
@show getFactor(retDFG, :x2x3x4f1).solverData.multihypo

# check for match
@test getFactor(fg, :x2x3x4f1).solverData.multihypo - getFactor(retDFG, :x2x3x4f1).solverData.multihypo |> norm < 1e-10
@test getFactor(fg, :x2x3x4f1).solverData.certainhypo - getFactor(retDFG, :x2x3x4f1).solverData.certainhypo |> norm < 1e-10

##
end


@testset "saving to and loading from FileDFG with nullhypo, eliminated, solveInProgress" begin
##

fg = generateGraph_Kaess()
getSolverParams(fg).attemptGradients = true

addVariable!(fg, :x4, ContinuousScalar)
addFactor!(fg, [:x2;:x3;:x4], LinearRelative(Normal()), multihypo=[1.0;0.6;0.4])
addFactor!(fg, [:x1;], Prior(Normal(10,1)), nullhypo=0.5)

solveTree!(fg)

#manually change a few fields to test if they are preserved
fa = getFactor(fg, :x2x3x4f1)
fa.solverData.eliminated = true
fa.solverData.solveInProgress = 1
fa.solverData.nullhypo = 0.5


saveFolder = "/tmp/dfg_test"
saveDFG(fg, saveFolder)

retDFG = initfg()
getSolverParams(retDFG).attemptGradients = true
loadDFG!(retDFG, saveFolder)
Base.rm(saveFolder*".tar.gz")

@test issetequal(ls(fg), ls(retDFG))
@test issetequal(lsf(fg), lsf(retDFG))

@show getFactor(fg, :x2x3x4f1).solverData.multihypo
@show getFactor(retDFG, :x2x3x4f1).solverData.multihypo

# check for match
@test isapprox(getFactor(fg, :x2x3x4f1).solverData.multihypo, getFactor(retDFG, :x2x3x4f1).solverData.multihypo)
@test isapprox(getFactor(fg, :x2x3x4f1).solverData.certainhypo, getFactor(retDFG, :x2x3x4f1).solverData.certainhypo)


fb = getFactor(retDFG, :x2x3x4f1)
@test fa == fb

fa = getFactor(fg, :x1f2)
fb = getFactor(retDFG, :x1f2)

@test fa == fb

##
end

