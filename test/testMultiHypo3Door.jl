# load requried packages
using IncrementalInference
using Test


## parameters

lm_prior_noise = 0.01
meas_noise = 0.25
odom_noise = 0.1
n_samples = 200

# initialize mean landmark locations
l0 = 0.0
l1 = 10.0
l2 = 20.0
l3 = 40.0

# "Ground-truth" robot poses
x0 = 0.0
x1 = 10.0
x2 = 20.0
x3 = 40.0

##

@testset "Basic 3 door, trinary multihypothesis tests..." begin

## Initialize empty factor graph
fg = initfg()
getSolverParams(fg).N = n_samples
getSolverParams(fg).gibbsIters = 5

# Place strong prior on locations of three "doors"
addVariable!(fg, :l0, ContinuousScalar, N=n_samples)
addFactor!(fg, [:l0], Prior(Normal(l0, lm_prior_noise)))

addVariable!(fg, :l1, ContinuousScalar, N=n_samples)
addFactor!(fg, [:l1], Prior(Normal(l1, lm_prior_noise)))

addVariable!(fg, :l2, ContinuousScalar, N=n_samples)
addFactor!(fg, [:l2], Prior(Normal(l2, lm_prior_noise)))

addVariable!(fg, :l3, ContinuousScalar, N=n_samples)
addFactor!(fg, [:l3], Prior(Normal(l3, lm_prior_noise)))


# Add first pose
addVariable!(fg, :x0, ContinuousScalar, N=n_samples)

# Make first "door" measurement
addFactor!(fg, [:x0; :l0; :l1; :l2; :l3], LinearRelative(Normal(0, meas_noise)), multihypo=[1.0; (1/4 for _=1:4)...])

# Add second pose
addVariable!(fg, :x1, ContinuousScalar, N=n_samples)

# Gaussian transition model
addFactor!(fg, [:x0; :x1], LinearRelative(Normal(x1-x0, odom_noise)))

# Make second "door" measurement
addFactor!(fg, [:x1; :l0; :l1; :l2; :l3], LinearRelative(Normal(0, meas_noise)), multihypo=[1.0; (1/4 for _=1:4)...])

##

solveGraph!(fg)

##

# check there is enough likelihood in the right places
@test 0.1 < getBelief(fg, :x0)([l0])[1]
@test 0.1 < getBelief(fg, :x0)([l1])[1]
@test getBelief(fg, :x0)([l2])[1] < 0.3

@test getBelief(fg, :x1)([l0])[1] < 0.3
@test 0.1 < getBelief(fg, :x1)([l1])[1]
@test 0.1 < getBelief(fg, :x1)([l2])[1]


##

for i in 1:1
  solveGraph!(fg);
end

##


# check there is enough likelihood in the right places
@test 0.1 < getBelief(fg, :x0)([l0])[1]
@test 0.1 < getBelief(fg, :x0)([l1])[1]
@test getBelief(fg, :x0)([l2])[1] < 0.03

@test getBelief(fg, :x1)([l0])[1] < 0.03
@test 0.1 < getBelief(fg, :x1)([l1])[1]
@test 0.1 < getBelief(fg, :x1)([l2])[1]


## Add one more pose/odometry to invoke issue #236

# Add third pose
addVariable!(fg, :x2, ContinuousScalar, N=n_samples)
addFactor!(fg, [:x1; :x2], LinearRelative(Normal(x2-x1, odom_noise)))

## Solve graph

tree = solveTree!(fg)
# drawGraph(fg)
# drawTree(tree, show=true)

##

# check there is enough likelihood in the right places
@test 0.1 < getBelief(fg, :x0)([l0])[1]
@test 0.1 < getBelief(fg, :x0)([l1])[1]
@test getBelief(fg, :x0)([l2])[1] < 0.03

@test getBelief(fg, :x1)([l0])[1] < 0.03
@test 0.1 < getBelief(fg, :x1)([l1])[1]
@test 0.1 < getBelief(fg, :x1)([l2])[1]

dx = x2-x1
@test getBelief(fg, :x2)([l0 + dx])[1] < 0.03
@test 0.1 < getBelief(fg, :x2)([l1 + dx])[1]
@test 0.1 < getBelief(fg, :x2)([l2 + dx])[1]

##

# Add third pose
addVariable!(fg, :x3, ContinuousScalar, N=n_samples)
addFactor!(fg, [:x2; :x3], LinearRelative(Normal(x3-x2, odom_noise)))


# Make third "door" measurement
addFactor!(fg, [:x3; :l0; :l1; :l2; :l3], LinearRelative(Normal(0, meas_noise)), multihypo=[1.0; (1/4 for _=1:4)...])

##

solveGraph!(fg)

##

@test isapprox(mean(getBelief(fg, :x0))[1], x0; atol = 2.0)
@test isapprox(mean(getBelief(fg, :x1))[1], x1; atol = 2.0)
@test isapprox(mean(getBelief(fg, :x2))[1], x2; atol = 2.0)
@test isapprox(mean(getBelief(fg, :x3))[1], x3; atol = 2.0)

@test isapprox(mean(getBelief(fg, :l0))[1], l0; atol = 3.0)
@test isapprox(mean(getBelief(fg, :l1))[1], l1; atol = 3.0)
@test isapprox(mean(getBelief(fg, :l2))[1], l2; atol = 3.0)
@test isapprox(mean(getBelief(fg, :l3))[1], l3; atol = 3.0)

##

@error "diabling final tests for now, see #1570"

if false
# check the PPEs are the same
@test isapprox(getPPE(fg, :x0).suggested[1], x0; atol = 2.0)
@test isapprox(getPPE(fg, :x1).suggested[1], x1; atol = 2.0)
@test isapprox(getPPE(fg, :x2).suggested[1], x2; atol = 2.0)
@test isapprox(getPPE(fg, :x3).suggested[1], x3; atol = 2.0)

@test isapprox(getPPE(fg, :l0).suggested[1], l0; atol = 3.0)
@test isapprox(getPPE(fg, :l1).suggested[1], l1; atol = 3.0)
@test isapprox(getPPE(fg, :l2).suggested[1], l2; atol = 3.0)
@test isapprox(getPPE(fg, :l3).suggested[1], l3; atol = 3.0)
end


##


end



# using RoMEPlotting
# Gadfly.set_default_plot_size(35cm,25cm)
#
# plotBelief(fg, sortDFG(ls(fg, r"x")))


#
