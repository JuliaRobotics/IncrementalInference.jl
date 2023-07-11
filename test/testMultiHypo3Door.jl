# load requried packages
# using Revise
using IncrementalInference
using Test


## during dev its clear functionality is working with 8/10 quality (Test API makes it difficult to write deterministic only tests for 8/10 quality.)

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
# forcefully work with ccw.varValsAll to check the pointers are pointing to getVal.(variables)
getSolverParams(fg).graphinit = false

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
f1 = addFactor!(fg, [:x0; :l0; :l1; :l2; :l3], LinearRelative(Normal(0, meas_noise)), multihypo=[1.0; (1/4 for _=1:4)...])

# check pointers before init
a,b = IIF._checkVarValPointers(fg, getLabel(f1))
for i in 1:length(a)
  @test a[i] == b[i]
end

# do init (and check that the var pointers did not change)
doautoinit!(fg ,:l0)
doautoinit!(fg ,:l1)
doautoinit!(fg ,:l2)
doautoinit!(fg ,:l3)


# make sure approxConv is as expected
@test isInitialized.(fg, [:l0;:l1;:l2;:l3]) |> all
x0_beforeConv = getVal(fg, :x0) |> deepcopy

# do the computation
X0 = approxConvBelief(fg, getLabel(f1), :x0)
# smpls = sampleFactor(fg, f1.label,10)

# check that the x0 variable memory has not be changed
@test all(norm.(x0_beforeConv - getVal(fg, :x0)) .< 1e-10)

# specifically after approxConv to :x0
a_,b_ = IIF._checkVarValPointers(fg, getLabel(f1))
# deep copy on destination memory for x0, given just did approxConv to x0s
@test a_[1] != b_[1]
@test a_[2] == b_[2]
@test a_[3] == b_[3]
@test a_[4] == b_[4]
@test a_[5] == b_[5]


##

# should have four equal sized peaks at landmark locations
@test 0.1 < X0([l0])[1]
@test 0.1 < X0([l1])[1]
@test 0.1 < X0([l2])[1]
@test 0.1 < X0([l3])[1]


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

for i in 1:2
  solveGraph!(fg);
end

##


# check there is enough likelihood in the right places
@test 0.1 < getBelief(fg, :x0)([l0])[1]
@test 0.1 < getBelief(fg, :x0)([l1])[1]
@test getBelief(fg, :x0)([l2])[1] < 0.03

# @test getBelief(fg, :x1)([l0])[1] < 0.03 # why add this?
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
@test 0.05 < getBelief(fg, :x0)([l0])[1]
@test 0.05 < getBelief(fg, :x0)([l1])[1]

@test 0.05 < getBelief(fg, :x1)([l1])[1]
@test 0.05 < getBelief(fg, :x1)([l2])[1]

dx = x2-x1
@test 0.05 < getBelief(fg, :x2)([l1 + dx])[1]
@test 0.05 < getBelief(fg, :x2)([l2 + dx])[1]

if false
  @test getBelief(fg, :x0)([l2])[1] < 0.03
  @test getBelief(fg, :x1)([l0])[1] < 0.03
  @test getBelief(fg, :x2)([l0 + dx])[1] < 0.03
else
  @error("Suppressing parts of multihypo tests (stochastic pass or fail results in poor test quality")
end

##

# Add third pose
addVariable!(fg, :x3, ContinuousScalar, N=n_samples)
addFactor!(fg, [:x2; :x3], LinearRelative(Normal(x3-x2, odom_noise)))


# Make third "door" measurement
addFactor!(fg, [:x3; :l0; :l1; :l2; :l3], LinearRelative(Normal(0, meas_noise)), multihypo=[1.0; (1/4 for _=1:4)...])

##

solveGraph!(fg)

##

@error "must restore a few multimodal tests"
if false
@test isapprox(mean(getBelief(fg, :x0))[1], x0; atol = 3.0)
@test isapprox(mean(getBelief(fg, :x1))[1], x1; atol = 3.0)
# @test isapprox(mean(getBelief(fg, :x2))[1], x2; atol = 3.0)
@test isapprox(mean(getBelief(fg, :x3))[1], x3; atol = 3.0)

@test isapprox(mean(getBelief(fg, :l0))[1], l0; atol = 3.0)
@test isapprox(mean(getBelief(fg, :l1))[1], l1; atol = 3.0)
@test isapprox(mean(getBelief(fg, :l2))[1], l2; atol = 3.0)
@test isapprox(mean(getBelief(fg, :l3))[1], l3; atol = 3.0)

##


# check the PPEs are the same
@test isapprox(getPPE(fg, :x0).suggested[1], x0; atol = 3.0)
@test isapprox(getPPE(fg, :x1).suggested[1], x1; atol = 3.0)
@test isapprox(getPPE(fg, :x2).suggested[1], x2; atol = 3.0)
@test isapprox(getPPE(fg, :x3).suggested[1], x3; atol = 3.0)

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
