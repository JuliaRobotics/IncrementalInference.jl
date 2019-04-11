## load the libraries
# option load before IIF
# using Cairo, Fontfconfig, Gadfly

# the multimodal iSAM library
using IncrementalInference

# build some factor graph
fg = initfg()
addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0], Prior(Normal(0,1)))
addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0, :x1], LinearConditional(Normal(10.0,1)))
addVariable!(fg, :x2, ContinuousScalar)
mmo = MixtureLinearConditional([Rayleigh(3); Uniform(30,55)], Categorical([0.4; 0.6]))
addFactor!(fg, [:x1, :x2], mmo)


# show the factor graph
writeGraphPdf(fg, show=true)
# show the tree
tree = wipeBuildNewTree!(fg, drawpdf=true, show=true)


# solve the factor graph and show solving progress on tree in src/JunctionTree.jl
tree = batchSolve!(fg, drawpdf=true, show=true)


## building a new tree -- as per IIF.prepBatchTree(...)

const IIF = IncrementalInference

IIF.resetFactorGraphNewTree!(fg)

# Look at variable ordering used to build the Bayes net/tree
p = IIF.getEliminationOrder(fg, ordering=:qr)


fge = deepcopy(fg)

# Building Bayes net.
IIF.buildBayesNet!(fge, p)

# prep and build tree
tree = emptyBayesTree()
IIF.buildTree!(tree, fge, p)

# Find potential functions for each clique
cliq = tree.cliques[1] # start at the root
IIF.buildCliquePotentials(fg, tree, cliq);

IIF.drawTree(tree, show=true)

# println("Bayes Net")
# sleep(0.1)
#fid = open("bn.dot","w+")
#write(fid,to_dot(fge.bn))
#close(fid)


## can also show the Clique Association matrix by first importing Cairo, Fontconfig, Gadfly

cliq = tree.cliques[1]
cliq = whichCliq(tree, :x0) # where is :x0 a frontal variable
spyCliqMat()

tree = drawTree(tree, show=true, imgs=true)
