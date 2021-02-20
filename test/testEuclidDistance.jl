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
## Test zero with y-axis
points = [[0.0;100.0],[100.0;0.0]]
fg = IIF.generateCanonicalFG_EuclidDistance(points)
# solveTree!(fg)

# fg = initfg()
# # getSolverParams(fg).inflation=250.0

# addVariable!(fg, :x0, ContinuousEuclid{2}, N=N)
# addFactor!(fg, [:x0], Prior(MvNormal([100.0;0], [1;1.0])), graphinit=false)

# addVariable!(fg, :x1, ContinuousEuclid{2}, N=N)
# addFactor!(fg, [:x1], Prior(MvNormal([0.0;100.0], [1;1.0])), graphinit=false)

# addVariable!(fg, :l1, ContinuousEuclid{2}, N=N)
# addFactor!(fg, [:x0;:l1], EuclidDistance(Normal(100.0, 1.0)), graphinit=false)
# addFactor!(fg, [:x1;:l1], EuclidDistance(Normal(100.0, 1.0)), graphinit=false)

## 

eo = [:x2; :x1; :l1]

tree = buildTreeReset!(fg, eo)

# drawTree(tree, show=true)
# drawGraph(fg, show=true)


## what would clique solution produce as up message

# solveTree!(fg)

# @error "continue test dev with #1168"
#solve the clique in isolation
hist = solveCliq!(fg, tree, :x2; recordcliq=true)
printCliqHistorySummary(hist)

# the belief that would have been sent by this clique:
belief = IIF.getMessageBuffer(hist[11].csmc.cliq).upTx
# L1 = getCliqueData(getClique(tree, :x1)).messages.upTx.belief[:l1] |> manikde!

fnc_, csmc_ = repeatCSMStep!(hist, 5);
sfg = csmc_.cliqSubFg

plotKDE(sfg, :l1)

##

IIF._getCCW(sfg, :x2l1f1, :l1)

pts = approxConv(sfg, :x2l1f1, :l1 , skipSolve=false)
initManual!(sfg, :l1, pts)
pts = approxConv(sfg, :x1l1f1, :l1)

plotKDE(manikde!(pts, ContinuousEuclid{2}))

##

end




@testset "Euclid Distance Tests" begin
# using Random
# Random.seed!(84)
# N=100
##
points = [[100.0],]
fg = IIF.generateCanonicalFG_EuclidDistance(points)
solveTree!(fg)

@test isapprox(getPPE(fg, :x1).suggested[1], 100, atol=1)

pts = getBelief(fg, :l1) |> getPoints
N = size(pts, 2)

# TODO add similar tests to the rest
@test_broken 0.3*N < sum(isapprox.(pts,  0, atol=5)) < 0.7*N
@test_broken 0.3*N < sum(isapprox.(pts,200, atol=5)) < 0.7*N

# Do it manually with a big inflation
# IIF._getCCW(fg, :x1l1f1).inflation = 100.0 # never gets there
# IIF._getCCW(fg, :x1l1f1).inflation = 150.0 # few iters gets there
IIF._getCCW(fg, :x1l1f1).inflation = 200.0 # One almost, second good
pts = approxConv(fg, :x1l1f1, :l1)
initManual!(fg, :l1, pts)
# plotKDE(fg, ls(fg))

pts = approxConv(fg, :x1l1f1, :l1)
initManual!(fg, :l1, pts)
# plotKDE(fg, ls(fg))

@test 0.3*N < sum(isapprox.(pts,  0, atol=5)) < 0.7*N
@test 0.3*N < sum(isapprox.(pts,200, atol=5)) < 0.7*N


## Test zero with x-axis
points = [[100.0;0.0],]
fg = IIF.generateCanonicalFG_EuclidDistance(points)
solveTree!(fg)

## Test zero with y-axis
points = [[0.0;100.0],]
fg = IIF.generateCanonicalFG_EuclidDistance(points)
solveTree!(fg)

## Test zero with xy-axis 2 points
points = [[0.0;100.0],[100.0;0.0]]
fg = IIF.generateCanonicalFG_EuclidDistance(points)
solveTree!(fg)

## Test offsett with xy-axis 2 points
points = [[50.0;100.0],[100.0;50.0]]
fg = IIF.generateCanonicalFG_EuclidDistance(points; dist=50.0)
solveTree!(fg)
# plotKDE(fg, ls(fg))

## Manual init
points = [[0.0;100.0],[100.0;0.0]]
fg = IIF.generateCanonicalFG_EuclidDistance(points)
getSolverParams(fg).inflation=3.0

initManual!(fg, :x1, rand(MvNormal([100.,0], [1.,1]),N))
initManual!(fg, :x2, rand(MvNormal([0.,100], [1.,1]),N))

# init = MixtureModel([MvNormal([100.,100], [10.,10]),
#                        MvNormal([0.,0], [10.,10])],
#                        [0.5, 0.5])
init = MvNormal([25.,25], [1.,1])
initManual!(fg, :l1, rand(init,N))

# plotKDE(fg, ls(fg))

# normal 2 clique eliminationOrder 
eliminationOrder = [:l1; :x2; :x1]
# one clique eliminationOrder
eliminationOrder = [:l1; :x2; :x1]
tree,_ = solveTree!(fg; eliminationOrder)

end

## SolverPlotter debug
# Random.seed!(84)
# empty!(IIF.g_u0)
# empty!(IIF.g_r)
# Plots.scatter([getindex.(IIF.g_u0,1), getindex.(IIF.g_u0,2)], legend=nothing)
# Plots.scatter!([getindex.(IIF.g_r,1),getindex.(IIF.g_r,2)], legend=nothing)
# Plots.scatter([getindex.(IIF.g_r,1),getindex.(IIF.g_r,2)], legend=nothing)
# x = reshape(getindex.(IIF.g_r,1),100,:)
# y = reshape(getindex.(IIF.g_r,2),100,:)
# Plots.scatter(x[:,1:2:end],y[:,1:2:end], legend=nothing)
# Plots.scatter(x[:,1:2:end],y[:,1:2:end], legend=nothing)


## what would clique solution produce as up message
