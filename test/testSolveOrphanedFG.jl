# test forest of graphs can solve with CSM, specifically #518

using IncrementalInference
using Statistics
using Test



@testset "Test forest of orphaned graphs" begin

fg = initfg()
addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(0,0.1)))
addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0;:x1], LinearConditional(Normal(10,0.1)))
addVariable!(fg, :x2, ContinuousScalar)
addFactor!(fg, [:x1;:x2], LinearConditional(Normal(10,0.1)))

addVariable!(fg, :x10, ContinuousScalar)
addFactor!(fg, [:x10;], Prior(Normal()))
addVariable!(fg, :x11, ContinuousScalar)
addFactor!(fg, [:x10;:x11], LinearConditional(Normal(-10,1.0)))
addVariable!(fg, :x12, ContinuousScalar)
addFactor!(fg, [:x11;:x12], LinearConditional(Normal(-10,1.0)))

# dfgplot(fg)
# getSolverParams(fg).drawtree = true
# getSolverParams(fg).showtree = true
# solve factor graph with two orphaned components
vo = Symbol[:x12, :x2, :x0, :x11, :x1, :x10]
tree, smt, hist = solveTree!(fg, variableOrder=vo)

# test tree will have two different root nodes
@test getVariableOrder(tree) == vo

@test getParent(tree, getCliq(tree, :x1)) |> length == 0
@test getParent(tree, getCliq(tree, :x10)) |> length == 0

@test getChildren(tree, getCliq(tree, :x1)) |> length == 1
@test getChildren(tree, getCliq(tree, :x10)) |> length == 1

@test getChildren(tree, getCliq(tree, :x2)) |> length == 0
@test getChildren(tree, getCliq(tree, :x12)) |> length == 0


## Test the numerical values are correct

@test getKDE(fg, :x0) |> getPoints  |> mean |> abs < 1.0
@test (getKDE(fg, :x1) |> getPoints  |> mean) - 10 |> abs < 2.0
@test (getKDE(fg, :x2) |> getPoints  |> mean) - 20 |> abs < 3.0

@test getKDE(fg, :x10) |> getPoints  |> mean |> abs < 2.0
@test (getKDE(fg, :x11) |> getPoints  |> mean) + 10 |> abs < 4.0
@test (getKDE(fg, :x12) |> getPoints  |> mean) + 20 |> abs < 5.0



# using RoMEPlotting
# Gadfly.set_default_plot_size(35cm, 25cm)
#
# plotKDE(fg, ls(fg))


end
