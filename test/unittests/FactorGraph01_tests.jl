## Initial setup
# using Revise
using Juno
using IncrementalInference
using RoME
using DistributedFactorGraphs
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

# FYI deep copies not required, since no per factor specific data in pp
addFactor!(dfg, [:x1, :x2], pp, autoinit=false)
addFactor!(dfg, [:x2, :x3], pp, autoinit=false)
addFactor!(dfg, [:x1, :l1], deepcopy(p2br), autoinit=false)
# TODO: Ask why it can't resolve here.
addFactor!(dfg, [:x2, :l1], deepcopy(p2br), autoinit=false)
addFactor!(dfg, [:x3, :l2], deepcopy(p2br), autoinit=false)

# Show it
DFG.toDotFile(dfg, "/tmp/testRmMarg.dot")

# Test with a copy
dfgPrime = deepcopy(dfg)
elimOrder = getEliminationOrder(dfgPrime)
# Forcing elimination order to match the documentation
elimOrder = [:l1, :l2, :x1, :x2, :x3]

# tree for elim order: [:l1, :l2, :x1, :x2, :x3]
# :l1 <-- :x1, :x2
# :x1 <-- :x2
# :x2 <-- :x3
# :l2 <-- :x3
# :x3 <-- /-/

buildBayesNet!(dfgPrime, elimOrder)
# Reproduce example
# IncrementalInference.rmVarFromMarg(dfgPrime, lm[1], )
DFG.toDotFile(dfgPrime, "/tmp/testRmMarg.dot")
# Assert that everything was eliminated and every variable has a BayesNetVertID
@test all(map(v -> solverData(v.dfgNode).eliminated, values(dfgPrime.g.vertices)))
# @test all(map(v -> getData(v).BayesNetVertID != nothing, DFG.getVariables(dfgPrime)))
# Assert that we have the expected Bayes tree
expectedBayesOutVertDict = Dict{Symbol, Vector{Symbol}}(:x2 => [:x3], :l1 => [:x1, :x2], :x3 => [], :l2 => [:x3], :x1 => [:x2])
for (vId,linkIds) in expectedBayesOutVertDict
    v = DFG.getVariable(dfgPrime, vId)
    @test setdiff(solverData(v).BayesNetOutVertIDs, linkIds) == []
end

# Now build the tree
# global tree = emptyBayesTree()
# buildTree!(tree, dfgPrime, p)

## End-to-end
@info "Time to build the Bayes Tree..."
tree = wipeBuildNewTree!(dfg, drawpdf=true, show=false)
drawTree(tree, show=true)


# TODO: Test
# Upward solve steps with clique state machine
# @info "Complete upward solve..."
# smtasks, ch = initInferTreeUp!(dfg, tree, drawtree=true, recordcliqs=[:x1;:x2;:x3] )




## Upward solve steps with clique state machine
# @info "Complete upward solve..."
# smtasks, ch = initInferTreeUp!(fgl, bt, N=N, drawtree=drawpdf, recordcliqs=recordcliqs, limititers=limititers, skipcliqids=skipcliqids )

# @info "solve leaf clique with single state machine"
# alltasks[i] = @async tryCliqStateMachineSolve!(fgl, treel, i, cliqHistories, drawtree=drawtree, N=N, limititers=limititers, recordcliqs=recordcliqs)
# manually

resetTreeCliquesForUpSolve!(tree)
setTreeCliquesMarginalized!(dfg, tree)



cliq = whichCliq(tree, :x1)



# history = cliqInitSolveUpByStateMachine!(dfg, tree, cliq, drawtree=true,
#                                          limititers=50, recordhistory=true )
# 0


children = Graphs.ExVertex[]
for ch in Graphs.out_neighbors(cliq, tree.bt)
  push!(children, ch)
end
prnt = getParent(tree, cliq)
csmc = CliqStateMachineContainer(dfg, initfg(), tree, cliq, prnt, children, false, true, true, getSolverParams(dfg))
statemachine = StateMachine{CliqStateMachineContainer}(next=testCliqCanRecycled_StateMachine) # isCliqUpSolved_StateMachine)


limititers = 100
recordhistory = true


statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory)
statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory)
statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory)
statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory)
statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory)
statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory)

statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory)

lsf(csmc.cliqSubFg, :x2)

writeGraphPdf(csmc.cliqSubFg, show=true, engine="neato")


statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory)
statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory)



csmc.cliqSubFg.nodeCounter



writeGraphPdf(csmc.cliqSubFg, show=true, engine="neato")




addFactor!(dfg, [:x1, :x2], pp, autoinit=false)
DFG.lsf(dfg, :x1x2f1)


0



















statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory)
statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory)
statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory)











DFG.getFactorIds(dfg)




cliqHistories = Dict{Int,Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}()
isTreeSolved(treel, skipinitialized=true)



# TODO: Nary tests for factors....



DFG.getFactor(dfg, :x1x2f1)



## See picture of upward clique association matrix

using Gadfly
cliq = whichCliq(tree, :x1)
spyCliqMat(cliq)



###
