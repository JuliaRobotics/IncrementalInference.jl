
using Test
using IncrementalInference


##

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


##

# deleteVariableSolverData!.(fg, ls(fg))
# @test listSolveKeys(fg) |> length == 0

##

getSolverParams(fg).dbg=true
# error("remember not 15")
getSolverParams(fg).limititers=30

# tree = buildTreeReset!(fg)


# ##

# stuff = solveCliqUp!(fg, tree, 2, solveKey=:testSolveKey)
# stuff = solveCliqUp!(fg, tree, 3, solveKey=:testSolveKey)
# stuff = solveCliqUp!(fg, tree, 1, solveKey=:testSolveKey)


##

smtasks = Task[]
solveTree!(fg, smtasks=smtasks, verbose=true, recordcliqs=ls(fg) )# , solveKey=:testSolveKey )
hists = fetchCliqHistoryAll!(smtasks)

##

using RoMEPlotting
Gadfly.set_default_plot_size(35cm,25cm)

##

plotSLAM2D(fg, solveKey=:testSolveKey)

##

end

#