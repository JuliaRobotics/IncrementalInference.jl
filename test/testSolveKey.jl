
using Test
using IncrementalInference


##

@testset "Basic test of belief prediction on alternate solveKey" begin

##


fg = initfg()

addVariable!(fg, :a, ContinuousScalar)
addVariable!(fg, :b, ContinuousScalar)

addFactor!(fg, [:a], Prior(Normal(10,1)), graphinit=false)
addFactor!(fg, [:a;:b], LinearRelative(Normal(10,1)), graphinit=false)


deleteVariableSolverData!(fg, :a, :default)
deleteVariableSolverData!(fg, :b, :default)

##

pts = sampleFactor(fg, :af1, 100)

IIF.setDefaultNodeData!(getVariable(fg, :a), 0, 100, 1, solveKey=:testSolveKey, 
                        initialized=false, varType=ContinuousScalar(), dontmargin=false)
#

initManual!(fg, :a, pts, :testSolveKey)

@test !isInitialized(fg, :a)
@test isInitialized(fg, :a, :testSolveKey)

##

IIF.setDefaultNodeData!(getVariable(fg, :b), 0, 100, 1, solveKey=:testSolveKey, 
                        initialized=false, varType=ContinuousScalar(), dontmargin=false)
#


@test !(:default in listSolveKeys(getVariable(fg, :a)))
@test !(:default in listSolveKeys(getVariable(fg, :b)))

@test (:testSolveKey in listSolveKeys(getVariable(fg, :a)))
@test (:testSolveKey in listSolveKeys(getVariable(fg, :b)))


##

doautoinit!(fg, :b, solveKey=:testSolveKey)


##


@test !(:default in listSolveKeys(getVariable(fg, :a)))
@test !(:default in listSolveKeys(getVariable(fg, :b)))

@test (:testSolveKey in listSolveKeys(getVariable(fg, :a)))
@test (:testSolveKey in listSolveKeys(getVariable(fg, :b)))


##

@test isapprox( calcPPE(fg, :a, solveKey=:testSolveKey).suggested[1], 10, atol=1)
@test isapprox( calcPPE(fg, :b, solveKey=:testSolveKey).suggested[1], 20, atol=1)


##

end


@testset "test solve with unique solveKey, see #1219" begin

##

fg = initfg()
getSolverParams(fg).graphinit=false

addVariable!(fg, :a, ContinuousScalar)
addVariable!(fg, :b, ContinuousScalar)
addVariable!(fg, :c, ContinuousScalar)
addVariable!(fg, :d, ContinuousScalar)
addVariable!(fg, :e, ContinuousScalar)

addFactor!(fg, [:a], Prior(Normal()))
addFactor!(fg, [:a;:b], LinearRelative(Normal(10, 1)))
addFactor!(fg, [:b;:c], LinearRelative(Normal(10, 1)))
addFactor!(fg, [:c;:d], LinearRelative(Normal(10, 1)))
addFactor!(fg, [:d;:e], LinearRelative(Normal(10, 1)))

getSolverParams(fg).graphinit=true


##

# getSolverParams(fg).limititers=30
solveTree!(fg, solveKey=:testSolveKey )

##

@test isapprox( calcPPE(fg, :a, solveKey=:testSolveKey).suggested[1], 0, atol=2)
@test isapprox( calcPPE(fg, :b, solveKey=:testSolveKey).suggested[1], 10, atol=2)
@test isapprox( calcPPE(fg, :c, solveKey=:testSolveKey).suggested[1], 20, atol=2)
@test isapprox( calcPPE(fg, :d, solveKey=:testSolveKey).suggested[1], 30, atol=2)
@test isapprox( calcPPE(fg, :e, solveKey=:testSolveKey).suggested[1], 40, atol=2)


##

# using RoMEPlotting
# Gadfly.set_default_plot_size(35cm,25cm)

# ##

# plotKDE(fg, ls(fg), solveKey=:testSolveKey)

##

end

#