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



@testset "test upward clique message range density behavior" begin

##

N=100
fg = initfg()
# getSolverParams(fg).inflation=250.0

addVariable!(fg, :x0, ContinuousEuclid{2}, N=N)
addFactor!(fg, [:x0], Prior(MvNormal([100.0;0], [1;1.0])), graphinit=false)

addVariable!(fg, :x1, ContinuousEuclid{2}, N=N)
addFactor!(fg, [:x1], Prior(MvNormal([0.0;100.0], [1;1.0])), graphinit=false)

addVariable!(fg, :l1, ContinuousEuclid{2}, N=N)
addFactor!(fg, [:x0;:l1], EuclidDistance(Normal(100.0, 1.0)), graphinit=false)
addFactor!(fg, [:x1;:l1], EuclidDistance(Normal(100.0, 1.0)), graphinit=false)

## 

eo = [:x1; :x0; :l1]

tree = buildTreeReset!(fg, eo)

## what would clique solution produce as up message

# solveTree!(fg)

# @error "continue test dev with #1168"
#solve the clique in isolation
stuff = solveCliq!(fg, tree, :x1; recordcliq=true)

# the belief that would have been sent by this clique:
belief = IIF.getMessageBuffer(stuff[11].csmc.cliq).upTx

## still need to make sure numerical results are fine..., first must resolve #1168

# L1 = getCliqueData(getClique(tree, :x1)).messages.upTx.belief[:l1] |> manikde!


##

end



#