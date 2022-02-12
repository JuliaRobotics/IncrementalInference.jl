# test forest of graphs can solve with CSM, specifically #518

using IncrementalInference
using Statistics
using Test
using TensorCast


@testset "Test forest of orphaned graphs" begin

fg = initfg()
addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(0,0.1)))
addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0;:x1], LinearRelative(Normal(10,0.1)))
addVariable!(fg, :x2, ContinuousScalar)
addFactor!(fg, [:x1;:x2], LinearRelative(Normal(10,0.1)))

addVariable!(fg, :x10, ContinuousScalar)
addFactor!(fg, [:x10;], Prior(Normal()))
addVariable!(fg, :x11, ContinuousScalar)
addFactor!(fg, [:x10;:x11], LinearRelative(Normal(-10,1.0)))
addVariable!(fg, :x12, ContinuousScalar)
addFactor!(fg, [:x11;:x12], LinearRelative(Normal(-10,1.0)))

# plotDFG(fg)
# getSolverParams(fg).drawtree = true
# getSolverParams(fg).showtree = true
# solve factor graph with two orphaned components
vo = Symbol[:x12, :x2, :x0, :x11, :x1, :x10]
tree = solveTree!(fg, eliminationOrder=vo)

# test tree will have two different root nodes
@test getEliminationOrder(tree) == vo

@test getParent(tree, getClique(tree, :x1)) |> length == 0
@test getParent(tree, getClique(tree, :x10)) |> length == 0

@test getChildren(tree, getClique(tree, :x1)) |> length == 1
@test getChildren(tree, getClique(tree, :x10)) |> length == 1

@test getChildren(tree, getClique(tree, :x2)) |> length == 0
@test getChildren(tree, getClique(tree, :x12)) |> length == 0


## Test the numerical values are correct

pts_ = getBelief(fg, :x0) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test  pts  |> mean |> abs < 1.0
pts_ = getBelief(fg, :x1) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test (pts  |> mean) - 10 |> abs < 2.0
pts_ = getBelief(fg, :x2) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test (pts  |> mean) - 20 |> abs < 3.0

pts_ = getBelief(fg, :x10) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test  pts  |> mean |> abs < 2.0
pts_ = getBelief(fg, :x11) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test (pts  |> mean) + 10 |> abs < 4.0
pts_ = getBelief(fg, :x12) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test (pts  |> mean) + 20 |> abs < 5.0



# using RoMEPlotting
# Gadfly.set_default_plot_size(35cm, 25cm)
#
# plotKDE(fg, ls(fg))


end
