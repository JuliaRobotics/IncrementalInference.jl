using Revise
using IncrementalInference
using RoME
using DistributedFactorGraphs
import DistributedFactorGraphs.GraphsJl
const DFGGraphs = DistributedFactorGraphs.GraphsJl
using Test


dfg = initfg()
addVariable!(dfg, :x1, Pose2)
addVariable!(dfg, :x2, Pose2)
addVariable!(dfg, :x3, Pose2)
addVariable!(dfg, :l1, Pose2)
addVariable!(dfg, :l2, Pose2)

prior = PriorPose2( MvNormal([10; 10; pi/6.0], Matrix(Diagonal([0.1;0.1;0.05].^2))))
addFactor!(dfg, [:x1], prior )

pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))

addFactor!(dfg, [:x1, :x2], deepcopy(pp))
addFactor!(dfg, [:x2, :x3], deepcopy(pp))
addFactor!(dfg, [:x1, :l1], deepcopy(p2br))
# TODO: Ask why it can't resolve here.
addFactor!(dfg, [:x2, :l1], deepcopy(p2br))
addFactor!(dfg, [:x3, :l2], deepcopy(p2br))

# Show it
DFGGraphs.toDotFile(dfg, "/tmp/testRmMarg.dot")

# Test with a copy
dfgPrime = deepcopy(dfg)
elimOrder = IncrementalInference.getEliminationOrder(dfgPrime)
IncrementalInference.buildBayesNet!(dfgPrime, elimOrder)
# Reproduce example
# IncrementalInference.rmVarFromMarg(dfgPrime, lm[1], )
DFGGraphs.toDotFile(dfgPrime, "/tmp/testRmMarg.dot")
