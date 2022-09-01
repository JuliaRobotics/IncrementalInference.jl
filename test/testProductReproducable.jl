# test for conv and product repeatability


using Test
using Statistics
using IncrementalInference

##

@testset "forward backward convolutions and products sequence" begin

fg = initfg()

addVariable!(fg, :a, ContinuousScalar)
addVariable!(fg, :b, ContinuousScalar)
addVariable!(fg, :c, ContinuousScalar)
addVariable!(fg, :d, ContinuousScalar)
addVariable!(fg, :e, ContinuousScalar)

addFactor!(fg, [:a], Prior(Normal()))
addFactor!(fg, [:a;:b], LinearRelative(Normal(10, 1)))
addFactor!(fg, [:b;:c], LinearRelative(Normal(10, 1)))
addFactor!(fg, [:c;:d], LinearRelative(Normal(10, 1)))
addFactor!(fg, [:d;:e], LinearRelative(Normal(10, 1)))

initAll!(fg)

tree = solveTree!(fg)


@test (Statistics.mean(getPoints(getBelief(fg, :a)))- 0 |> abs) < 3
@test (Statistics.mean(getPoints(getBelief(fg, :b)))-10 |> abs) < 4
@test (Statistics.mean(getPoints(getBelief(fg, :c)))-20 |> abs) < 4
@test (Statistics.mean(getPoints(getBelief(fg, :d)))-30 |> abs) < 5
@test (Statistics.mean(getPoints(getBelief(fg, :e)))-40 |> abs) < 5

@test 0.3 < Statistics.std(getPoints(getBelief(fg, :a))) < 2
@test 0.5 < Statistics.std(getPoints(getBelief(fg, :b))) < 4
@test 0.9 < Statistics.std(getPoints(getBelief(fg, :c))) < 6
@test 1.2 < Statistics.std(getPoints(getBelief(fg, :d))) < 7
@test 1.5 < Statistics.std(getPoints(getBelief(fg, :e))) < 8


# drawTree(tree, show=true)
# using RoMEPlotting
# plotKDE(fg, ls(fg))
# spyCliqMat(tree, :b)

end


@testset "Basic back and forth convolution over LinearRelative should spread" begin

fg = initfg()

addVariable!(fg, :a, ContinuousScalar)
addVariable!(fg, :b, ContinuousScalar)

addFactor!(fg, [:a;:b], LinearRelative(Normal(10, 1)), graphinit=false)

initVariable!(fg, :a, randn(1,100))
initVariable!(fg, :b, 10 .+randn(1,100))

A = getBelief(fg, :a)
B = getBelief(fg, :b)
# plotKDE(fg, [:a; :b])

# repeat many times to ensure the means stay put and covariances spread out
for i in 1:10
  pts = approxConv(fg, :abf1, :b)
  B_ = manikde!(ContinuousScalar, pts)
  # plotKDE([B_; B])
  initVariable!(fg, :b, B_)

  pts = approxConv(fg, :abf1, :a)
  A_ = manikde!(ContinuousScalar, pts)
  # plotKDE([A_; A])
  initVariable!(fg, :a, A_)
end

A_ = getBelief(fg, :a)
B_ = getBelief(fg, :b)
# plotKDE([A_; B_; A; B])

@test (Statistics.mean(getPoints(A)) |> abs) < 1
@test (Statistics.mean(getPoints(A_))|> abs) < 2

@test (Statistics.mean(getPoints(B)) -10 |> abs) < 1
@test (Statistics.mean(getPoints(B_))-10 |> abs) < 2

@test Statistics.std(getPoints(A)) < 2
@test 3 < Statistics.std(getPoints(A_))

@test Statistics.std(getPoints(B)) < 2
@test 3 < Statistics.std(getPoints(B_))

##

end



##
