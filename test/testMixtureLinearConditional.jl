using IncrementalInference
using Test

#
@testset "test simple MixtureLinearConditional" begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0], Prior(Normal(0.0,0.1)))

mlr = MixtureRelative(LinearRelative(I), [Normal(-1.0, 0.1), Normal(1.0, 0.1)], Categorical([0.5; 0.5]))

# test serialization while we are here
pmlr = convert(PackedMixtureRelative, mlr)
umlr = convert(MixtureRelative, pmlr)

@test mlr.mechanics == umlr.mechanics
@test mlr.components == umlr.components
@test mlr.diversity == umlr.diversity

mlr = MixtureLinearConditional([Normal(-1.0, 0.1), Normal(1.0, 0.1)], Categorical([0.5; 0.5]))

addFactor!(fg, [:x0,:x1], mlr)

# To look at your factor graph
# if false
# using GraphPlot
# using DistributedFactorGraphs
# dfgplot(fg)
# end

tree, smt, hist = solveTree!(fg)

btd = getBelief(getVariable(fg, :x0))
@test isapprox(mean(getKDEfit(btd,distribution=Normal)), 0.0; atol=0.1) 

@test isapprox(std(getKDEfit(btd,distribution=Normal)), 0.1; atol=0.05) 

btd = getBelief(getVariable(fg, :x1))
pts = getPoints(btd)
pts_p = pts[pts .>= 0]
pts_n = pts[pts .< 0]

nfit_p = fit(Normal, pts_p)
@test isapprox(mean(nfit_p), 1.0; atol=0.1)
@test isapprox(std(nfit_p), 0.14; atol=0.05) #TODO confirm the correct value and tolerance

nfit_n = fit(Normal, pts_n)
@test isapprox(mean(nfit_n), -1.0; atol=0.1)
@test isapprox(std(nfit_n), 0.14; atol=0.05) #TODO confirm the correct value and tolerance

# To look at your results
# if false
#   using RoMEPlotting
#   plotKDE(fg, ls(fg))
# end

end
