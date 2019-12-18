# IIF #485 --

using Test
using Logging
using DistributedFactorGraphs
using IncrementalInference



logger = SimpleLogger(stdout, Logging.Debug)
global_logger(logger)
dfg = initfg() #LightDFG{SolverParams}(solverParams=SolverParams())
# Add some nodes.
v1 = addVariable!(dfg, :a, ContinuousScalar, labels = [:POSE])
v2 = addVariable!(dfg, :b, ContinuousScalar, labels = [:POSE])
v3 = addVariable!(dfg, :c, ContinuousScalar, labels = [:LANDMARK])
f1 = addFactor!(dfg, [:a; :b], LinearConditional(Normal(50.0,2.0)) )
f2 = addFactor!(dfg, [:b; :c], LinearConditional(Normal(50.0,2.0)) )

drawGraph(dfg, show=true)

getSolverParams(dfg).limititers = 50
getSolverParams(dfg).dbg = true
getSolverParams(dfg).dbg = true
getSolverParams(dfg).drawtree = true
getSolverParams(dfg).showtree = true

# tree = wipeBuildNewTree!(dfg)
#
# using Gadfly, Cairo, Fontconfig
#
# drawTree(tree, show=true, imgs=true)

tree, smtasks, hist = solveTree!(dfg, recordcliqs=ls(dfg))




#
