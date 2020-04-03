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
l2 = 40.0

# "Ground-truth" robot poses
x0 = 0.0
x1 = 10.0
x2 = 20.0
x3 = 40.0


@testset "2door basic binary multihypothesis tests..." begin

## Initialize empty factor graph
fg = initfg()
fg.solverParams.N = n_samples
# Place strong prior on locations of three "doors"
addVariable!(fg, :l0, ContinuousScalar, N=n_samples)
addFactor!(fg, [:l0], Prior(Normal(l0, lm_prior_noise)))

addVariable!(fg, :l1, ContinuousScalar, N=n_samples)
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

## Solve graph

tree, smt, hist = solveTree!(fg)
# drawGraph(fg)
# drawTree(tree, show=true)


@test abs(getKDEMean(getKDE(fg, :x0))[1]-x0) < 2.0
@test abs(getKDEMean(getKDE(fg, :x1))[1]-x1) < 2.0
@test abs(getKDEMean(getKDE(fg, :x2))[1]-x2) < 2.0

@test abs(getKDEMean(getKDE(fg, :l0))[1]-l0) < 3.0
@test abs(getKDEMean(getKDE(fg, :l1))[1]-l1) < 3.0


0
end

#
# using RoMEPlotting
# Gadfly.set_default_plot_size(35cm,25cm)
#
# plotKDE(fg, [:l0;:l1])
# plotKDE(fg, [:x0;:x1;:x2])


@testset "3door basic trinary multihypo test..." begin

# Initialize empty factor graph
fg = initfg()

# Place strong prior on locations of three "doors"
addVariable!(fg, :l0, ContinuousScalar, N=n_samples)
addFactor!(fg, [:l0], RoME.Prior(Normal(l0, lm_prior_noise)))

addVariable!(fg, :l1, ContinuousScalar, N=n_samples)
addFactor!(fg, [:l1], RoME.Prior(Normal(l1, lm_prior_noise)))

addVariable!(fg, :l2, ContinuousScalar, N=n_samples)
addFactor!(fg, [:l2], RoME.Prior(Normal(l2, lm_prior_noise)))

# Add first pose
addVariable!(fg, :x0, ContinuousScalar, N=n_samples)

# Make first "door" measurement
addFactor!(fg, [:x0; :l0; :l1; :l2], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/3.0; 1.0/3.0; 1.0/3.0])

# Add second pose
addVariable!(fg, :x1, ContinuousScalar, N=n_samples)

# # Gaussian transition model
addFactor!(fg, [:x0; :x1], LinearConditional(Normal(x1-x0, odom_noise)))

# # Make second "door" measurement
addFactor!(fg, [:x1; :l0; :l1; :l2], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/3.0; 1.0/3.0; 1.0/3.0])

# # Add third pose
addVariable!(fg, :x2, ContinuousScalar, N=n_samples)

# # Gaussian transition model
addFactor!(fg, [:x1; :x2], LinearConditional(Normal(x2-x1, odom_noise)))

# Add fourth pose
addVariable!(fg, :x3, ContinuousScalar, N=n_samples)

# Add odometry transition and new landmark sighting
addFactor!(fg, [:x2, :x3], LinearConditional(Normal(x3-x2, odom_noise)))
addFactor!(fg, [:x3; :l0; :l1; :l2], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/3.0; 1.0/3.0; 1.0/3.0])


ensureAllInitialized!(fg)

drawGraph(fg)
solveTree!(fg)



@test abs(getKDEMean(getKDE(fg, :x0))[1]-x0) < 2.0
@test abs(getKDEMean(getKDE(fg, :x1))[1]-x1) < 2.0
@test abs(getKDEMean(getKDE(fg, :x2))[1]-x2) < 2.0

@test abs(getKDEMean(getKDE(fg, :l0))[1]-l0) < 3.0
@test abs(getKDEMean(getKDE(fg, :l1))[1]-l1) < 3.0
@test abs(getKDEMean(getKDE(fg, :l2))[1]-l2) < 3.0


end


# using RoMEPlotting
# Gadfly.set_default_plot_size(35cm,25cm)
#
# plotKDE(fg, [:l0;:l1;:l2])
# plotKDE(fg, [:x0;:x1;:x2])


#
