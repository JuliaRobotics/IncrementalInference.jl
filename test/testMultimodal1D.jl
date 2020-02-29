
using DistributedFactorGraphs
using IncrementalInference

using Test


##==============================================================================
## Setup
##==============================================================================

n_samples = 100
graphinit = false

# noises
lm_prior_noise = 1
meas_noise = 1
odom_noise = 1

# initialize mean landmark locations
l1 = -30.0
l2 = 30.0
l3 = -40.0

p_meas = 0.5
p_map = 0.5

graphinit = false

##==============================================================================
# Initialize empty factor graph

fg = initfg()
fg.solverParams.N = n_samples

# lp landmark prior information
# lm landmark measurement

addVariable!(fg, :lp1, ContinuousScalar, autoinit=graphinit, N=n_samples)
addFactor!(fg, [:lp1], Prior(Normal(l1, lm_prior_noise)))

addVariable!(fg, :lp2, ContinuousScalar, autoinit=graphinit, N=n_samples)
addFactor!(fg, [:lp2], Prior(Normal(l2, lm_prior_noise)))

addVariable!(fg, :x1, ContinuousScalar, autoinit=graphinit, N=n_samples)

addVariable!(fg, :lm2, ContinuousScalar, autoinit=graphinit, N=n_samples)
addFactor!(fg, [:x1; :lm2; :lp2], LinearConditional(Normal(20., meas_noise)), multihypo=[1.0; p_meas; p_map])


addVariable!(fg, :lm1, ContinuousScalar, autoinit=graphinit, N=n_samples)
addFactor!(fg, [:x1; :lm1; :lp1], LinearConditional(Normal(-20., meas_noise)), multihypo=[1.0; p_meas; p_map])

#weak lp1 lm1 relation
addFactor!(fg, [:lp1; :lm1], LinearConditional(Normal(0., 100.)))



using RoMEPlotting
using Plots


#
