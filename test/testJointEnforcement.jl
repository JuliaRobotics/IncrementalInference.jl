# test case to enforce consistency in joint gibbs

using Test
using IncrementalInference


@testset "test case with disjoint clique joint subgraph" begin

## test case with disjoint clique joint subgraph

fg = initfg()

addVariable!(fg, :x0, ContinuousEuclid{2})
addVariable!(fg, :x1, ContinuousEuclid{2})
addVariable!(fg, :x2, ContinuousEuclid{2})

initManual!(fg, :x0, randn(2,100))
initManual!(fg, :x1, randn(2,100) .+ 10)
initManual!(fg, :x2, randn(2,100) .+ 20)

addFactor!(fg , [:x0; :x1], LinearRelative(MvNormal([10.0;10], diagm([1.0;1]))))
addFactor!(fg , [:x1; :x2], LinearRelative(MvNormal([10.0;10], diagm([1.0;1]))))

# setPPE!(fg, :x2)
# fg[:x2]


##

addVariable!(fg, :x3, ContinuousEuclid{2})
addFactor!( fg, [:x2; :x3], EuclidDistance(Normal(10, 1)) )

##

addFactor!( fg, [:x0; :x3], EuclidDistance(Normal(30, 1)), graphinit=false )


##

# drawGraph(fg, show=true)


## test shortest path

# first path
tp1_ = [ :x0;:x0x1f1;:x1;:x1x2f1;:x2]
tp2_ = [ :x0;:x0x3f1;:x3;:x2x3f1;:x2]

pth = findShortestPathDijkstra(fg, :x0, :x2)
@test pth == tp1_ || pth == tp2_

pth = findShortestPathDijkstra(fg, :x0, :x2, typeFactors=[LinearRelative;])
@test pth == tp1_

# different path
pth = findShortestPathDijkstra(fg, :x0, :x2, typeFactors=[EuclidDistance;])
@test pth == tp2_


## use a specific solve order

vo = [:x3; :x1; :x2; :x0] # getEliminationOrder(fg)
tree, _, = solveTree!(fg, variableOrder=vo);

##

# drawTree(tree, show=true)

## get up message from child clique

msgBuf = IIF.getMessageBuffer(getClique(tree, :x3))
msg = msgBuf.upTx


## Child clique subgraph

# each clique has a subgraph
cliq2 = getClique(tree,:x3)
cfg2 = buildCliqSubgraph(fg, cliq2)

drawGraph(cfg2, show=true)

##


separators = getCliqSeparatorVarIds(cliq2)

allClasses = IIF._findSubgraphsFactorType( cfg2, msg.diffJoints, separators )

upmsgpriors = _generateSubgraphMsgPriors( cfg2, msg, allClasses)


##





##

prs = _generateSubgraphMsgPriors(cfg2, msg, allClasses)


##


mb = IIF.getMessageBuffer(getClique(tree, :x0))

mb.upRx[2].diffJoints[1].variables
mb.upRx[2].diffJoints[1].likelihood


##

@show findShortestPathDijkstra(fg, :x0,:x2)

@show findShortestPathDijkstra(fg, :x1,:x3)



##

isPathFactorsHomogeneous(fg, :x0, :x2)


##

fg[:x4]
fg[:x2]

##


# vo = getEliminationOrder(fg)

tree = resetBuildTree!(fg)
# tree = resetBuildTreeFromOrder!(fg, [:x0;:x1;:x3;:x2])

drawTree(tree, show=true)



##





## check which path between separators has homogeneous factors


isHom, ftyps = isPathFactorsHomogeneous(fg, :x0, :x2)


_sft = selectFactorType(fg, :x0, :x2) 
sft = _sft()

typeof(sft).name == ftyps[1]

getindex(Main, ftyps[1])


##  dev


##


retlist = addLikelihoodsDifferentialCHILD!([:x2; :x0], cfg2)


retlist[1][2]




##

drawGraph(tfg, show=true)


##

getManifolds(LinearRelative)





#
