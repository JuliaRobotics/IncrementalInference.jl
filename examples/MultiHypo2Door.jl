# load requried packages
using RoME
using Test


## parameters

lm_prior_noise = 0.01
meas_noise = 0.25
odom_noise = 0.1
n_samples = 100

# initialize mean landmark locations
l0 = 0.0
l1 = 10.0
l2 = 40.0

# "Ground-truth" robot poses
x0 = 0.0
x1 = 10.0
x2 = 20.0
x3 = 40.0

## Initialize empty factor graph
fg = initfg()

# Place strong prior on locations of three "doors"
addVariable!(fg, Symbol("l0"), ContinuousScalar, N=n_samples)
addFactor!(fg, [:l0], Prior(Normal(l0, lm_prior_noise)))

addVariable!(fg, Symbol("l1"), ContinuousScalar, N=n_samples)
addFactor!(fg, [:l1], Prior(Normal(l1, lm_prior_noise)))


# Add first pose
addVariable!(fg, :x0, ContinuousScalar, N=n_samples)

# Make first "door" measurement
# addFactor!(fg, [:x0; :l0], LinearConditional(Normal(0, meas_noise)))
addFactor!(fg, [:x0; :l0; :l1], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/2.0; 1.0/2.0])


# Add second pose
addVariable!(fg, :x1, ContinuousScalar, N=n_samples)

# Gaussian transition model
addFactor!(fg, [:x0; :x1], LinearConditional(Normal(x1-x0, odom_noise)))

# Make second "door" measurement
# addFactor!(fg, [:x1; :l1], LinearConditional(Normal(0, meas_noise)) )
addFactor!(fg, [:x1; :l0; :l1], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/2.0; 1.0/2.0])




## Add one more pose/odometry to invoke issue #236

# Add third pose
addVariable!(fg, :x2, ContinuousScalar, N=n_samples)
addFactor!(fg, [:x1; :x2], LinearConditional(Normal(x2-x1, odom_noise)))


# Add fourth pose
# addVariable!(fg, :x3, ContinuousScalar, N=n_samples)

# Add odometry transition and new landmark sighting
# addFactor!(fg, [:x2, :x3], LinearConditional(Normal(2, odom_noise)))
# addFactor!(fg, [:x3; :l0; :l1; :l2], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/3.0; 1.0/3.0; 1.0/3.0])

## Do some debugging
ensureAllInitialized!(fg)

##

writeGraphPdf(fg, show=true)

tree = wipeBuildNewTree!(fg, drawpdf=true, show=true)

## Solve graph
tree = batchSolve!(fg)

# tree = wipeBuildNewTree!(fg, drawpdf=true, show=true)



## Plotting functions below

using RoMEPlotting


pl = plotKDE(fg, [:x0;:x1])

pl = plotKDE(fg, [:x0;:x1;:x2])
pl |> PNG("/tmp/test.png")

pl = plotKDE(fg, [:l0; :l1])

spyCliqMat(tree, :l0)

spyCliqMat(tree, :x2)



## specialized debugging


stuff = treeProductUp(fg, tree, :l0, :x0)
plotKDE(manikde!(stuff[1], (:Euclid,)) )


## Do one clique inference only

tree = wipeBuildNewTree!(fg, drawpdf=true, show=true)
urt = doCliqInferenceUp!(fg, tree, :l0, false, iters=1, drawpdf=true)
upmsgs = urt.keepupmsgs

plotKDE([upmsgs[:x0]; upmsgs[:l1]; upmsgs[:x1]], c=["red";"green";"blue"])


## swap iteration order

getData(tree.cliques[2]).itervarIDs = [5;7;3;1]

inferOverTree!(fg, tree)

## manually build the iteration scheme for second clique

# iter order:  x0, x1, l0, l1

stuff = treeProductUp(fg, tree, :l0, :x0)
X0 = manikde!(stuff[1], (:Euclid,))
plotKDE([X0; getKDE(fg, :x0)], c=["red";"green"])
setValKDE!(fg, :x0, X0)


stuff = treeProductUp(fg, tree, :l0, :x1)
X1 = manikde!(stuff[1], (:Euclid,))
plotKDE([X1; getKDE(fg, :x1)], c=["red";"green"])
setValKDE!(fg, :x1, X1)


stuff = treeProductUp(fg, tree, :l0, :l0)
L0 = manikde!(stuff[1], (:Euclid,))
plotKDE([L0; getKDE(fg, :l0)], c=["red";"green"])
setValKDE!(fg, :l0, L0)


stuff = treeProductUp(fg, tree, :l0, :l1)
L1 = manikde!(stuff[1], (:Euclid,))
plotKDE([L1; getKDE(fg, :l1)], c=["red";"green"])
setValKDE!(fg, :l1, L1)



## Reconstruct individual steps for broader clique factor selection

# Cliq 2:
# L0,X0,L1,X1
# x ,  ,  ,
# x ,x ,x ,
# x ,  ,x ,x
#   ,x ,  ,x
#   ,  ,x ,

# choose iteration order (priors last): :x0, :x1, :l0, :l1

# for new initialization format
# cliq 2: init :l0, :l1 directly from priors, them proceed with regular order
# cliq 1: initialize :x2 from incoming message singleton and proceed with regular order


# get factors for :x0 in clique2:
# :x0l0l1f1, :x0x1f1
ptsX0, = predictbelief(fg, :x0, [:x0l0l1f1; :x0x1f1])
X0 = manikde!(ptsX0, (:Euclid,))
plotKDE([X0; getKDE(fg, :x0)], c=["red";"green"])
setValKDE!(fg, :x0, X0)

# get factors for :x1 in clique2:
# :x1l0l1f1, :x0x1f1
ptsX1, = predictbelief(fg, :x1, [:x1l0l1f1, :x0x1f1])
X1 = manikde!(ptsX1, (:Euclid,))
plotKDE([X1; getKDE(fg, :x1)], c=["red";"green"])
setValKDE!(fg, :x1, X1)


# get factors for :l0
# :x0l0l1f1, :x1l0l1f1, :l0f1
ptsL0, = predictbelief(fg, :l0, [:x0l0l1f1, :x1l0l1f1, :l0f1])
L0 = manikde!(ptsL0, (:Euclid,))
plotKDE([L0; getKDE(fg, :l0)], c=["red";"green"])
setValKDE!(fg, :l0, L0)


# get factors for :l1
# :x0l0l1f1, :x1l0l1f1, :l1f1
ptsL1, = predictbelief(fg, :l1, [:x0l0l1f1, :x1l0l1f1, :l1f1])
L1 = manikde!(ptsL1, (:Euclid,))
plotKDE([L1; getKDE(fg, :l1)], c=["red";"green"])
setValKDE!(fg, :l1, L1)



## double check the upmessages are stored properly


##

plotLocalProduct(fg, :x0)
plotLocalProduct(fg, :x1)
plotLocalProduct(fg, :l0)
plotLocalProduct(fg, :l1)


##


ensureAllInitialized!(fg)


tree = wipeBuildNewTree!(fg, drawpdf=true, show=true)
cliqorder = getCliqOrderUpSolve(tree)


spyCliqMat(cliqorder[end])



## Development zone

# treel = deepcopy(tree)
# fgl = deepcopy(fg)
# cliql = deepcopy(cliq)




##  develop better factor selection method

varlist = [:l0; :x0; :l1; :x1]
getFactorsAmongVariablesOnly(fg, varlist)




#
