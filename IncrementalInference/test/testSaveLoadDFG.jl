
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

@show DFG.getRecipehyper(fg, :x2x3x4f1).multihypo
@show DFG.getRecipehyper(retDFG, :x2x3x4f1).multihypo

# check for match
@test DFG.getRecipehyper(fg, :x2x3x4f1).multihypo - DFG.getRecipehyper(retDFG, :x2x3x4f1).multihypo |> norm < 1e-10

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
fa.state.eliminated = true
fa.hyper.nullhypo = 0.5


saveFolder = "/tmp/dfg_test"
saveDFG(fg, saveFolder)

retDFG = initfg()
getSolverParams(retDFG).attemptGradients = true
loadDFG!(retDFG, saveFolder)
Base.rm(saveFolder*".tar.gz")

@test issetequal(ls(fg), ls(retDFG))
@test issetequal(lsf(fg), lsf(retDFG))

@show DFG.getRecipehyper(fg, :x2x3x4f1).multihypo
@show DFG.getRecipehyper(retDFG, :x2x3x4f1).multihypo

# check for match
@test isapprox(DFG.getRecipehyper(fg, :x2x3x4f1).multihypo, DFG.getRecipehyper(retDFG, :x2x3x4f1).multihypo)


fb = getFactor(retDFG, :x2x3x4f1)
@test fa == fb

fa = getFactor(fg, :x1f2)
fb = getFactor(retDFG, :x1f2)

@test fa == fb

##
end

