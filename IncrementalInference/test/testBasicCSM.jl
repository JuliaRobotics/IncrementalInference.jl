# IIF #485 --

# using Revise

using Test
using Logging
using Statistics
using DistributedFactorGraphs
using IncrementalInference
using TensorCast


@testset "test basic three variable graph with prior" begin

VAR1 = :a
VAR2 = :b
VAR3 = :c

# logger = SimpleLogger(stdout, Logging.Debug)
# global_logger(logger)
dfg = initfg() #LocalDFG{SolverParams}(solverParams=SolverParams())
# Add some nodes.
v1 = addVariable!(dfg, VAR1, ContinuousScalar, tags = [:POSE])
v2 = addVariable!(dfg, VAR2, ContinuousScalar, tags = [:POSE])
v3 = addVariable!(dfg, VAR3, ContinuousScalar, tags = [:LANDMARK])
f1 = addFactor!(dfg, [VAR1; VAR2], LinearRelative(Normal(50.0,2.0)) )
f2 = addFactor!(dfg, [VAR2; VAR3], LinearRelative(Normal(50.0,2.0)) )

addFactor!(dfg, [VAR1], Prior(Normal()))

# drawGraph(dfg, show=true)


# tree = buildTreeReset!(dfg)
# # drawTree(tree, show=true)
#
# getCliqFactors(tree, VAR3)
# getCliqFactors(tree, VAR1)

initAll!(dfg)


# cliq= getClique(tree, VAR3)
# getCliqueData(cliq)
#
# cliq= getClique(tree, VAR1)
# getCliqueData(cliq)



getSolverParams(dfg).limititers = 50
# getSolverParams(dfg).drawtree = true
# getSolverParams(dfg).showtree = true
# getSolverParams(dfg).dbg = true
## getSolverParams(dfg).async = true


tree = solveTree!(dfg) #, recordcliqs=ls(dfg))

pts_ = getBelief(dfg, :c) |> getPoints
TensorCast.@cast pts[i,j] := pts_[j][i]

@test 70 < Statistics.mean(pts) < 130

# #
# using Gadfly, Cairo, Fontconfig
# drawTree(tree, show=true, imgs=true)

end




#
