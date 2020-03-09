# from #605

# using Revise

using IncrementalInference
using Test


@testset "test mixture prior" begin

# init graph
fg = initfg()
N = 100
fg.solverParams.N = N;

# add first variable
addVariable!(fg, :x0, ContinuousScalar)

# add bi-modal mixture prior
Prior0 = MixturePrior([Normal(-5.0,1.0); Normal(0.0,1.0)], Categorical([0.5;0.5]))
addFactor!(fg, [:x0], Prior0)

smpls, lb = getSample(Prior0, N)

# should be a balance of particles
@test sum(lb .== 1) - sum(lb .== 2) |> abs < 0.3*N
@test sum(smpls .< -2.5) - sum(-2.5 .< smpls) |> abs < 0.3*N

# solve
solveTree!(fg)

marginalPts = getKDE(fg, :x0) |> getPoints

# check solver solution consistent too
@test sum(marginalPts .< -2.5) - sum(-2.5 .< marginalPts) |> abs < 0.3*N


end


# using RoMEPlotting
#
# # plot the results
# plotKDE(fg, :x0)


#
