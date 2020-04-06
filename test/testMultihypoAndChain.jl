using Test
using IncrementalInference


@testset "test multihypo chain example (see #462)..." begin

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
addFactor!(fg, [:x1;:l1;:l1_0], LinearConditional(Normal(l1-x1, lnoise)), multihypo=[1;1/2;1/2])
addFactor!(fg, [:x2;:l2;:l2_0], LinearConditional(Normal(l2-x2, lnoise)), multihypo=[1;1/2;1/2])
addFactor!(fg, [:x1;:x2], LinearConditional(Normal(0, Onoise)))

tree, smt, hist = solveTree!(fg)

# drawTree(tree, show=true)

# expect x1 x2 to have at least one mode at 0

@test getPPE(fg, :x1).suggested[1] - x1 |> abs < 1
@test getPPE(fg, :x2).suggested[1] - x2 |> abs < 1

@test getPPE(fg, :l1).suggested[1] - l1 |> abs < 1
@test getPPE(fg, :l2).suggested[1] - l2 |> abs < 1

# l1_0, l2_0 should be nearby around l1 and l2, but cannot confirm 100%
@test getPPE(fg, :l1_0).suggested[1] - l1 |> abs < 10
@test getPPE(fg, :l2_0).suggested[1] - l2 |> abs < 10


end



# using RoMEPlotting
# Gadfly.set_default_plot_size(35cm,25cm)
#
# plotKDE(fg, [:l1;:l2])
# plotKDE(fg, [:l1_0;:l2_0])
# plotKDE(fg, [:x1;:x2])


#
