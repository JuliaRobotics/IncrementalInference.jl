# IIF #485 --

# using Revise

using Test
using Logging
using Statistics
using DistributedFactorGraphs
using IncrementalInference


@testset "test basic three variable graph with prior" begin

VAR1 = :a
VAR2 = :b
VAR3 = :c

logger = SimpleLogger(stdout, Logging.Debug)
global_logger(logger)
dfg = initfg() #LightDFG{SolverParams}(solverParams=SolverParams())
# Add some nodes.
v1 = addVariable!(dfg, VAR1, ContinuousScalar, labels = [:POSE])
v2 = addVariable!(dfg, VAR2, ContinuousScalar, labels = [:POSE])
v3 = addVariable!(dfg, VAR3, ContinuousScalar, labels = [:LANDMARK])
f1 = addFactor!(dfg, [VAR1; VAR2], LinearConditional(Normal(50.0,2.0)) )
f2 = addFactor!(dfg, [VAR2; VAR3], LinearConditional(Normal(50.0,2.0)) )

addFactor!(dfg, [VAR1], Prior(Normal()))

# drawGraph(dfg, show=true)


# tree = wipeBuildNewTree!(dfg)
# # drawTree(tree, show=true)
#
# getCliqFactors(tree, VAR3)
# getCliqFactors(tree, VAR1)

ensureAllInitialized!(dfg)


# cliq= getCliq(tree, VAR3)
# getData(cliq)
#
# cliq= getCliq(tree, VAR1)
# getData(cliq)



getSolverParams(dfg).limititers = 50
# getSolverParams(dfg).drawtree = true
# getSolverParams(dfg).showtree = true
# getSolverParams(dfg).dbg = true
## getSolverParams(dfg).async = true


tree, smtasks, hist = solveTree!(dfg) #, recordcliqs=ls(dfg))


@test 70 < Statistics.mean(getKDE(dfg, :c) |> getPoints) < 130

# #
# using Gadfly, Cairo, Fontconfig
# drawTree(tree, show=true, imgs=true)

end




#
