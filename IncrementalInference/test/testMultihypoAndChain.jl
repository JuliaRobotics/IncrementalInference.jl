using Test
using IncrementalInference
using Random

##

@testset "test basic multihypo" begin

## A simple multihypo example

Random.seed!(42) # The answer to reproducable noise

fg = LocalDFG(solverParams=SolverParams(graphinit=false, gibbsIters=5, spreadNH=5.0))

pRight = 0.99
pWrong = 0.01
pr_noise = 0.01
od_noise = 0.1
lm_noise = 0.01

# true positions
# x0 at 0
# x1 at 1
# l1 at 1
# l2 at 2

#x0 prior
addVariable!(fg, :x0, ContinuousScalar)
prpo = Normal(0.0, pr_noise)
addFactor!(fg, [:x0], Prior(Normal(rand(prpo), pr_noise)))

#l1 and l2
addVariable!(fg, :l1, ContinuousScalar, tags=[:LANDMARK])
addVariable!(fg, :l2, ContinuousScalar, tags=[:LANDMARK])

#x0 to l1 or l2
p2ln = Normal(1.0, lm_noise)
p2p = LinearRelative(Normal(rand(p2ln), lm_noise))
addFactor!(fg, [:x0; :l1; :l2], p2p, multihypo = [1, pRight, pWrong])
# addFactor!(fg, [:x0; :l1], p2p) #this one used for sanity check

#x0 to x1
addVariable!(fg, :x1, ContinuousScalar)
pp = Normal(1.0, od_noise)
addFactor!(fg, [:x0,:x1], LinearRelative(Normal(rand(pp), od_noise)))

#x1 to l1 or l2
p2ln = Normal(0.0, lm_noise)
p2p = LinearRelative(Normal(rand(p2ln), lm_noise))
addFactor!(fg, [:x1; :l1; :l2], p2p, multihypo = [1, pRight, pWrong])
# addFactor!(fg, [:x1; :l1], p2p) #this one used for sanity check

#x1 to l2 or l1
p2ln = Normal(1.0, lm_noise)
p2p = LinearRelative(Normal(rand(p2ln), lm_noise))
addFactor!(fg, [:x1; :l2; :l1], p2p, multihypo = [1, pRight, pWrong])
# addFactor!(fg, [:x1; :l2], p2p) #this one used for sanity check

##

# prescribe an elimination order to get a single clique
eo = [:l2,:x1,:x0,:l1]
# fg.solverParams.graphinit=true
smtasks = Task[]
tree = solveTree!(fg, eliminationOrder=eo) #, smtasks=smtasks, recordcliqs=ls(fg));


# hists = fetchCliqHistoryAll!(smtasks)

# plotKDE(fg, ls(fg))

##

@test isapprox(calcMeanMaxSuggested(fg, :x0, :default).suggested[], 0, atol = 0.2) 
@test isapprox(calcMeanMaxSuggested(fg, :x1, :default).suggested[], 1, atol = 0.2) 
@test isapprox(calcMeanMaxSuggested(fg, :l1, :default).suggested[], 1, atol = 0.2) 

L2 = getBelief(fg, :l2)
npts = length(getPoints(L2))
pts = [2.0.+0.1*randn(1) for _ in 1:npts]
L2_ = manikde!(ContinuousScalar, pts)

# test that there is at least a mode present
@test mmd(L2_, L2, ContinuousScalar) < 1e-3

##

end

@testset "test multihypo chain example (see #462)..." begin

##

l1 = -10.0
l2 = +10.0
lnoise = 1.0
x1 = 0
x2 = 0
Onoise = 0.1

fg = initfg()

addVariable!(fg, :x1, ContinuousScalar)
addVariable!(fg, :x2, ContinuousScalar)
addVariable!(fg, :l1, ContinuousScalar)
addVariable!(fg, :l1_0, ContinuousScalar)
addVariable!(fg, :l2, ContinuousScalar)
addVariable!(fg, :l2_0, ContinuousScalar)

# priors on two landmarks only
addFactor!(fg, [:l1], Prior(Normal(l1, lnoise)))
addFactor!(fg, [:l2], Prior(Normal(l2, lnoise)))

# relative constraints
addFactor!(fg, [:x1;:l1;:l1_0], LinearRelative(Normal(l1-x1, lnoise)), multihypo=[1;1/2;1/2])
addFactor!(fg, [:x2;:l2;:l2_0], LinearRelative(Normal(l2-x2, lnoise)), multihypo=[1;1/2;1/2])
addFactor!(fg, [:x1;:x2], LinearRelative(Normal(0, Onoise)))

tree = solveTree!(fg)

# drawTree(tree, show=true)

# expect x1 x2 to have at least one mode at 0

@test calcMeanMaxSuggested(fg, :x1, :default).suggested[1] - x1 |> abs < 1.2
@test calcMeanMaxSuggested(fg, :x2, :default).suggested[1] - x2 |> abs < 1.2

@test calcMeanMaxSuggested(fg, :l1, :default).suggested[1] - l1 |> abs < 1.2
@test calcMeanMaxSuggested(fg, :l2, :default).suggested[1] - l2 |> abs < 1.2

@test calcMeanMaxSuggested(fg, :l1_0, :default).suggested[1] - l1 |> abs < 10
@test calcMeanMaxSuggested(fg, :l2_0, :default).suggested[1] - l2 |> abs < 10

##

end



# using RoMEPlotting
# Gadfly.set_default_plot_size(35cm,25cm)
#
# plotKDE(fg, [:l1;:l2])
# plotKDE(fg, [:l1_0;:l2_0])
# plotKDE(fg, [:x1;:x2])


#
