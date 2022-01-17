# using Revise

using DistributedFactorGraphs
using IncrementalInference
using Test

##

@testset "test CSM runaway on upsolve, (issue 427)" begin

## parameters

lm_prior_noise = 0.1
meas_noise = 0.25
odom_noise = 0.1
n_samples = 100

# initialize mean landmark locations
l1 = 50.0
l2 = -50.0
l3 = 80.0

# "Ground-truth" robot poses
x1 = 0.0
x2 = 0.0
x3 = 0.0

## Initialize empty factor graph
fg = initfg()


addVariable!(fg, Symbol("l1"), ContinuousScalar, N=n_samples)
addFactor!(fg, [:l1], Prior(Normal(l1, lm_prior_noise)))

addVariable!(fg, Symbol("l2"), ContinuousScalar, N=n_samples)
addFactor!(fg, [:l2], Prior(Normal(l2, lm_prior_noise)))

addVariable!(fg, Symbol("l1_0"), ContinuousScalar, N=n_samples)
addVariable!(fg, Symbol("l2_0"), ContinuousScalar, N=n_samples)


# Add first pose
addVariable!(fg, :x1, ContinuousScalar, N=n_samples)

addFactor!(fg, [:x1; :l1; :l1_0], LinearRelative(Normal(40., meas_noise)), multihypo=[1.0; 1.0/2.0; 1.0/2.0])


# Add second pose
addVariable!(fg, :x2, ContinuousScalar, N=n_samples)

# Gaussian transition model
addFactor!(fg, [:x1; :x2], LinearRelative(Normal(0., odom_noise)))

# Make second "door" measurement
# addFactor!(fg, [:x1; :l1], LinearRelative(Normal(0, meas_noise)) )
addFactor!(fg, [:x2; :l2; :l2_0], LinearRelative(Normal(-40., meas_noise)), multihypo=[1.0; 1.0/2.0; 1.0/2.0])

# drawGraph(fg)

# initAll!(fg)

approxConv(fg, :x1l1l1_0f1, :l1_0)
approxConv(fg, :x1l1l1_0f1, :x1)
# doautoinit!(fg, :l1_0)
# doautoinit!(fg, :x1)


## Run solver
getSolverParams(fg).limititers = 30 # previous runaway CSM issue due to excessive limits on autoinit.
# getSolverParams(fg).dbg = false
# getSolverParams(fg).async = false
# getSolverParams(fg).drawtree = false
# getSolverParams(fg).showtree = false

tree = solveTree!(fg, recordcliqs=ls(fg))


# drawGraph(fg)
# fetchAssignTaskHistoryAll!(tree, smt)
# printCliqHistorySummary(tree, :l1)
# getTreeCliqsSolverHistories(fg, tree)
#
# csmAnimate(fg, tree, [:l1;:l2])
#
# # Base.rm("/tmp/caesar/csmCompound/out.ogv")
# run(`ffmpeg -r 10 -i /tmp/caesar/csmCompound/csm_%d.png -c:v libtheora -vf fps=25 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -q 10 /tmp/caesar/csmCompound/out.ogv`)
# @async run(`totem /tmp/caesar/csmCompound/out.ogv`)


# using RoMEPlotting
# plotKDE(fg, ls(fg))

##

end



#
