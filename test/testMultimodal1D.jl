
## Continue with loading packages

using DistributedFactorGraphs
using IncrementalInference
using TensorCast
using Test

##==============================================================================
## Setup
##==============================================================================

n_samples = 100
graphinit = false

# noises
lm_prior_noise = 1.0
meas_noise = 1.0
odom_noise = 1.0

# initialize mean landmark locations
l1 = -30.0
l2 = 30.0
# l3 = -40.0

p_meas = 0.4
p_map = 1 - p_meas

graphinit = false

##==============================================================================
# Initialize empty factor graph

@testset "test multihypo 1D..." begin

##

fg = initfg()
# fg.solverParams.graphinit = graphinit
fg.solverParams.N = n_samples
fg.solverParams.gibbsIters = 6
fg.solverParams.spreadNH = 0.3

# lp landmark prior information
# lm landmark measurement

addVariable!(fg, :lp1, ContinuousScalar)
addFactor!(fg, [:lp1], Prior(Normal(l1, lm_prior_noise)), graphinit=graphinit)

addVariable!(fg, :lp2, ContinuousScalar)
addFactor!(fg, [:lp2], Prior(Normal(l2, lm_prior_noise)), graphinit=graphinit)

addVariable!(fg, :x1, ContinuousScalar)

addVariable!(fg, :lm2, ContinuousScalar)
addFactor!(fg, [:x1; :lm2; :lp2], LinearRelative(Normal(20., meas_noise)), multihypo=[1.0; p_meas; p_map], graphinit=graphinit)

addVariable!(fg, :lm1, ContinuousScalar)
addFactor!(fg, [:x1; :lm1; :lp1], LinearRelative(Normal(-20., meas_noise)), multihypo=[1.0; p_meas; p_map], graphinit=graphinit)

#weak lp1 lm1 relation to nudge one of two symmetric options
# addFactor!(fg, [:lp1; :lm1], LinearRelative(Normal(0., 100.)))


# initAll!(fg)

# drawGraph(fg, show=true)
# getSolverParams(fg).drawtree = false
# getSolverParams(fg).showtree = false

# getSolverParams(fg).dbg = true
# getSolverParams(fg).multiproc = false

##

varor = [:x1, :lm1, :lm2, :lp1, :lp2]
# tree = buildTreeReset!(fg, varor)
tree = solveTree!(fg, eliminationOrder = varor);
# tree = solveTree!(fg, eliminationOrder = varor);

# plotKDE(fg, ls(fg))

##

# X should be at one of two modes
pts_ = getPoints(getBelief(fg, :x1))
@cast pts[i,j] := pts_[j][i]
@test 0.7*getSolverParams(fg).N < sum(-20 .< pts[:] .< 0) + sum(0 .< pts[:] .< 20)

pts_ = getPoints(getBelief(fg, :lp1))
@cast pts[i,j] := pts_[j][i]
@test 0.7*getSolverParams(fg).N < sum(-38 .< pts[:] .< -28)

pts_ = getPoints(getBelief(fg, :lp2))
@cast pts[i,j] := pts_[j][i]
@test 0.7*getSolverParams(fg).N < sum(28 .< pts[:] .< 38)

pts_ = getPoints(getBelief(fg, :lm1))
@cast pts[i,j] := pts_[j][i]
@test 0.1*getSolverParams(fg).N < sum(-38 .< pts[:] .< -25)

pts_ = getPoints(getBelief(fg, :lm2))
@cast pts[i,j] := pts_[j][i]
@test 0.1*getSolverParams(fg).N < sum(25 .< pts[:] .< 38)

# @test 0.2*getSolverParams(fg).N < sum(-18 .< getPoints(getBelief(fg, :lm1))[:] .< -5) ||
      # 0.2*getSolverParams(fg).N < sum(5 .< getPoints(getBelief(fg, :lm2))[:] .< 15)

##


end


## Debug plotting below

#
# using RoMEPlotting
# Gadfly.set_default_plot_size(35cm, 20cm)
#
#
# # tree = solveTree!(fg)
# varIds = [ :x1, :lp1, :lp2, :lm1, :lm2]
# pkde = plotKDE(fg, varIds)
#
# 0
#
# # # using Plots
# #
# # pmm = StatsPlots.plot(pkde.layers[2].mapping[:x],pkde.layers[2].mapping[:y], label = string(varIds[1]))
# # for i = 3:6
# #     StatsPlots.plot!(pmm, pkde.layers[i].mapping[:x],pkde.layers[i].mapping[:y], label=string(varIds[i-1]))
# # end
# # plot!(pmm, title = "MM-iSAM", xlims = (-50, 60), xticks = -50:10:50)
#
#
#
#
#
# p1 = plotLocalProduct(fg, :lm2)
# p2 = plotLocalProduct(fg, :lm1)
# p3 = plotLocalProduct(fg, :lp2)
# p4 = plotLocalProduct(fg, :lp1)
# p5 = plotLocalProduct(fg, :x1)
# h0 = hstack(p1,p2,p3,p4,p5)
#
# # treeProductUp(fg, tree, :x1, :x1)
#
#
#
# stuff = localProduct(fg,:lm2)
# initVariable!(fg,:lm2, stuff[1]); p1 = plotKDE(stuff[1], title="lm2")
#
# stuff = localProduct(fg,:lm1)
# initVariable!(fg,:lm1, stuff[1]); p2 = plotKDE(stuff[1], title="lm1")
#
# stuff = localProduct(fg,:lp2)
# initVariable!(fg,:lp2, stuff[1]); p3 = plotKDE(stuff[1], title="lp2")
#
# stuff = localProduct(fg,:lp1)
# initVariable!(fg,:lp1, stuff[1]); p4 = plotKDE(stuff[1], title="lp1")
#
# stuff = localProduct(fg,:x1)
# initVariable!(fg,:x1, stuff[1]); p5 = plotKDE(stuff[1], title="x1")
#
# h1 = hstack(p1,p2,p3,p4,p5)
#
# vstack(h0,h1,h2,h3,h4,h5) |> PDF("/tmp/test_new.pdf", 35cm, 40cm)
#
#
# fg1 = initfg()
# loadDFG("/tmp/fix/lm2_1.tar.gz", Main, fg1)
#
# fg2 = initfg()
# loadDFG("/tmp/fix/lm1_2.tar.gz", Main, fg2)
#
# fg3 = initfg()
# loadDFG("/tmp/fix/lp2_3.tar.gz", Main, fg3)
#
# fg4 = initfg()
# loadDFG("/tmp/fix/lp1_4.tar.gz", Main, fg4)
#
# fg5 = initfg()
# loadDFG("/tmp/fix/x1_5.tar.gz", Main, fg5)
#
#
# fg6 = initfg()
# loadDFG("/tmp/fix/lm2_6.tar.gz", Main, fg6)
#
# fg7 = initfg()
# loadDFG("/tmp/fix/lm1_7.tar.gz", Main, fg7)
#
# fg8 = initfg()
# loadDFG("/tmp/fix/lp2_8.tar.gz", Main, fg8)
#
# fg9 = initfg()
# loadDFG("/tmp/fix/lp1_9.tar.gz", Main, fg9)
#
# fg10 = initfg()
# loadDFG("/tmp/fix/x1_10.tar.gz", Main, fg10)
#
# fg11 = initfg()
# loadDFG("/tmp/fix/lm2_11.tar.gz", Main, fg11)
#
# fg12 = initfg()
# loadDFG("/tmp/fix/lm1_12.tar.gz", Main, fg12)
#
# fg13 = initfg()
# loadDFG("/tmp/fix/lp2_13.tar.gz", Main, fg13)
#
# fg14 = initfg()
# loadDFG("/tmp/fix/lp1_14.tar.gz", Main, fg14)
#
# fg15 = initfg()
# loadDFG("/tmp/fix/x1_15.tar.gz", Main, fg15)
#
#
# h1 = hstack(plotKDE(fg1,:lm2),plotKDE(fg2,:lm1),plotKDE(fg3,:lp2),plotKDE(fg4,:lp1),plotKDE(fg5,:x1))
#
# h2 = hstack(plotKDE(fg6,:lm2),plotKDE(fg7,:lm1),plotKDE(fg8,:lp2),plotKDE(fg9,:lp1),plotKDE(fg10,:x1))
#
# h3 = hstack(plotKDE(fg11,:lm2),plotKDE(fg12,:lm1),plotKDE(fg13,:lp2),plotKDE(fg14,:lp1),plotKDE(fg15,:x1))
#
#
# vstack(h1,h2,h3) |> PDF("/tmp/test_new_fail.pdf",35cm,30cm)
#
#




#
