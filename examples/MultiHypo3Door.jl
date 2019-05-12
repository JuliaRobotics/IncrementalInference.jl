# load requried packages
using RoME
using Test


## parameters

lm_prior_noise = 0.1
meas_noise = 0.15
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
addVariable!(fg, :l0, ContinuousScalar, N=n_samples)
addFactor!(fg, [:l0], Prior(Normal(l0, lm_prior_noise)))

addVariable!(fg, :l1, ContinuousScalar, N=n_samples)
addFactor!(fg, [:l1], Prior(Normal(l1, lm_prior_noise)))

addVariable!(fg, :l2, ContinuousScalar, N=n_samples)
addFactor!(fg, [:l2], Prior(Normal(l2, lm_prior_noise)))

# Add first pose
addVariable!(fg, :x0, ContinuousScalar, N=n_samples)

# Make first "door" measurement
addFactor!(fg, [:x0; :l0; :l1; :l2], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/3.0; 1.0/3.0; 1.0/3.0])

# Add second pose
addVariable!(fg, :x1, ContinuousScalar, N=n_samples)

# Gaussian transition model
addFactor!(fg, [:x0; :x1], LinearConditional(Normal(x1-x0, odom_noise)))

# Make second "door" measurement
addFactor!(fg, [:x1; :l0; :l1; :l2], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/3.0; 1.0/3.0; 1.0/3.0])



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

writeGraphPdf(fg, show=true)

tree = wipeBuildNewTree!(fg, drawpdf=true, show=true)


## Solve graph
tree = batchSolve!(fg)



## Plotting functions below

using RoMEPlotting


plotKDE(fg, [:x0])

plotKDE(fg, [:x0;:x1])

plotKDE(fg, [:x0;:x1;:x2])

plotKDE(fg, [:l0;:l1])

plotKDE(fg, [:l2])




spyCliqMat(tree, :l0)

spyCliqMat(tree, :x2)



## swap iteration order

getData(tree.cliques[2]).itervarIDs = [9;7;1;3;5]

inferOverTree!(fg, tree)





## Reconstruct individual steps for broader clique factor selection

# Cliq 2:
# L0,X0,L1,X1
# x ,  ,  ,
# x ,x ,x ,
# x ,  ,x ,x
#   ,x ,  ,x
#   ,  ,x ,

# choose iteration order (priors last): :x1, :x0, :l2, :l1, :l0

# get factors for :x1 in clique2:
# :x1l0l1l2f1, :x0x1f1
plotLocalProduct(fg, :x1, sidelength=20cm)
ptsX1, = predictbelief(fg, :x1, [:x1l0l1l2f1, :x0x1f1])
X1 = manikde!(ptsX1, (:Euclid,))
plotKDE([X1; getKDE(fg, :x1)], c=["red";"green"])
setValKDE!(fg, :x1, X1)


# get factors for :x0 in clique2:
# :x0l0l1l2f1, :x0x1f1
plotLocalProduct(fg, :x0, sidelength=20cm)
ptsX0, = predictbelief(fg, :x0, [:x0l0l1l2f1; :x0x1f1])
X0 = manikde!(ptsX0, (:Euclid,))
plotKDE([X0; getKDE(fg, :x0)], c=["red";"green"])
setValKDE!(fg, :x0, X0)


# get factors for :l2
# :x0l0l1l2f1, :x1l0l1l2f1, :l2f1
plotLocalProduct(fg, :l2, sidelength=20cm)
ptsL2, = predictbelief(fg, :l2, [:x0l0l1l2f1, :x1l0l1l2f1, :l2f1])
L2 = manikde!(ptsL2, (:Euclid,))
plotKDE([L2; getKDE(fg, :l2)], c=["red";"green"])
setValKDE!(fg, :l2, L2)


# get factors for :l1
# :x0l0l1l2f1, :x1l0l1l2f1, :l1f1
plotLocalProduct(fg, :l1, sidelength=20cm)
ptsL1, = predictbelief(fg, :l1, [:x0l0l1l2f1, :x1l0l1l2f1, :l1f1])
L1 = manikde!(ptsL1, (:Euclid,))
plotKDE([L1; getKDE(fg, :l1)], c=["red";"green"])
setValKDE!(fg, :l1, L1)


# get factors for :l0
# :x0l0l1l2f1, :x1l0l1l2f1, :l0f1
plotLocalProduct(fg, :l0, sidelength=20cm)
ptsL0, = predictbelief(fg, :l0, [:x0l0l1l2f1, :x1l0l1l2f1, :l0f1])
L0 = manikde!(ptsL0, (:Euclid,))
plotKDE([L0; getKDE(fg, :l0)], c=["red";"green"])
setValKDE!(fg, :l0, L0)



## Now do root clique manually too


ptsX2, = predictbelief(fg, :x2, [:x1x2f1;])
X2 = manikde!(ptsX2, (:Euclid,))
plotKDE([X2; getKDE(fg, :x2)], c=["red";"green"])
setValKDE!(fg, :x2, X2)


upmsgX1 = deepcopy(getKDE(fg, :x1))
ptsX1 = approxConv(fg, :x1x2f1, :x1)
pX1 = manikde!(ptsX1, (:Euclid,))
X1 = manifoldProduct([upmsgX1; pX1], (:Euclid,))
plotKDE([X1; getKDE(fg, :x1)], c=["red";"green"])
setValKDE!(fg, :x1, X1)





##  Complete downward pass


# :x0l0l1l2f1, :x1l0l1l2f1, :l1f1
ptsL1, = predictbelief(fg, :l1, [:x0l0l1l2f1, :x1l0l1l2f1, :l1f1])
L1 = manikde!(ptsL1, (:Euclid,))
plotKDE([L1; getKDE(fg, :l1)], c=["red";"green"])
setValKDE!(fg, :l1, L1)


## quick debug tests

##

pts = approxConv(fg, :x1l0l1l2f1, :l1)
plotKDE(manikde!(pts, (:Euclid,)))

##

plotKDE(fg, :x0)
pts = approxConv(fg, :x0l0l1l2f1, :l1)
plotKDE(manikde!(pts, (:Euclid,)))




## Debug draft-tree based autoinit:
# x0, l2, l1, l0, with one multihypo factor and 3 priors
# xxxx
#  x
#   x
#    x

fg = initfg()
addVariable!(fg, :l0, ContinuousScalar, N=n_samples)
addFactor!(fg, [:l0], Prior(Normal(l0, lm_prior_noise)), autoinit=false)
addVariable!(fg, :l1, ContinuousScalar, N=n_samples)
addFactor!(fg, [:l1], Prior(Normal(l1, lm_prior_noise)), autoinit=false)
addVariable!(fg, :l2, ContinuousScalar, N=n_samples)
addFactor!(fg, [:l2], Prior(Normal(l2, lm_prior_noise)), autoinit=false)
addVariable!(fg, :x0, ContinuousScalar, N=n_samples)
addFactor!(fg, [:x0; :l0; :l1; :l2], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/3.0; 1.0/3.0; 1.0/3.0], autoinit=false)
addVariable!(fg, :x1, ContinuousScalar, N=n_samples)
addFactor!(fg, [:x0; :x1], LinearConditional(Normal(x1-x0, odom_noise)),autoinit=false)
addFactor!(fg, [:x1; :l0; :l1; :l2], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/3.0; 1.0/3.0; 1.0/3.0],autoinit=false)
addVariable!(fg, :x2, ContinuousScalar, N=n_samples)
addFactor!(fg, [:x1; :x2], LinearConditional(Normal(x2-x1, odom_noise)),autoinit=false)



writeGraphPdf(fg, show=true)

tree = wipeBuildNewTree!(fg, drawpdf=true, show=true)

spyCliqMat(tree, :l0)
spyCliqMat(tree, :x2)

cliq = tree.cliques[2]

## quick debug



##

manualinit!(fg, :l2, [:l2f1])
manualinit!(fg, :l1, [:l1f1])
manualinit!(fg, :l0, [:l0f1])

# regular procedure
ptsX1, = predictbelief(fg, :x1, [:x1l0l1l2f1])
X1 = manikde!(ptsX1, (:Euclid,))
plotKDE(X1)
# just because the init flag has not been set yet
setValKDE!(fg, :x1, X1)

# regular procedure
ptsX0, = predictbelief(fg, :x0, [:x0l0l1l2f1])
X0 = manikde!(ptsX0, (:Euclid,))
plotKDE(X0)
# just because the init flag has not been set yet
setValKDE!(fg, :x0, X0)

# regular procedure
ptsX2, = predictbelief(fg, :x2, [:x1x2f1])
X2 = manikde!(ptsX2, (:Euclid,))
plotKDE(X2)
# just because the init flag has not been set yet
setValKDE!(fg, :x2, X2)







## other debugging


plotLocalProduct(fg, :x0)
plotLocalProduct(fg, :x1)
plotLocalProduct(fg, :l1, sidelength=15cm)

##

stuff = treeProductUp(fg, tree, :x2, :x2)


##

ensureAllInitialized!(fg)


tree = wipeBuildNewTree!(fg, drawpdf=true, show=true)
cliqorder = getCliqOrderUpSolve(tree)


spyCliqMat(cliqorder[end])



#
