# test EuclidDistance

using IncrementalInference
using Test

##

@testset "test EuclidDistance on 1 dim" begin

##

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0], Prior(Normal()))

addVariable!(fg, :x1, ContinuousScalar)

eud = EuclidDistance(Normal(10,1))
addFactor!(fg, [:x0;:x1], eud)

##

tree, _, = solveTree!(fg)

##

@test isapprox(getPPE(fg, :x0).suggested[1], 0, atol=1)

pts = getBelief(fg, :x1) |> getPoints
N = size(pts, 2)

@test 0.3*N < sum( 5 .< pts )
@test 0.3*N < sum( pts .< -5 )
@test sum( -5 .< pts .< 5 ) < 0.1*N

##

end


@testset "test EuclidDistance on 2 dim" begin

##

fg = initfg()

addVariable!(fg, :x0, ContinuousEuclid{2})
addFactor!(fg, [:x0], Prior(MvNormal(zeros(2),diagm([1;1.0]))))

addVariable!(fg, :x1, ContinuousEuclid{2})

eud = EuclidDistance(Normal(10,1))
addFactor!(fg, [:x0;:x1], eud)

tree, _, = solveTree!(fg)

@test isapprox(getPPE(fg, :x0).suggested[1], 0, atol=1)
@test isapprox(getPPE(fg, :x0).suggested[1], 0, atol=1)

pts = getBelief(fg, :x1) |> getPoints
N = size(pts, 2)

pts .^= 2
@test 0.5*N < sum( 7 .< sqrt.(sum(pts, dims=1)) .< 13 )

##

end