# test EuclidDistance

using IncrementalInference
using Test
using TensorCast

# using GLMakie
# function plotPoints(
#   hode::HomotopyDensity
# ) 
#   getDimension(hode) < 2 && error("plotPoints only implemented for 2 or more dimensions, got $(getDimension(hode))")
#   pts_ = getPoints(hode; permute=false)
#   @cast pts[i][j] := pts_[j][i]

#   scatter(pts...)
# end


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

tree = solveGraph!(fg)

##

@test isapprox(0, mean(getBelief(getState(fg, :x0, :default)))[1], atol=1)

pts_ = getBelief(fg, :x1) |> getPoints
@cast pts[i,j] := pts_[j][i]
N = size(pts, 2)

@test 0.3*N < sum( 5 .< pts )
@test 0.3*N < sum( pts .< -5 )
@test sum( -5 .< pts .< 5 ) < 0.1*N

##
end


@testset "2D range-only donut posterior test, utils test" begin
##

fg = initfg()

addVariable!(fg, :x0, ContinuousEuclid{2})
addFactor!(fg, [:x0], Prior(MvNormal(zeros(2),diagm([1;1.0]))))

# sanity check so far
doautoinit!(fg, :x0)
X0_ = getBelief(fg, :x0, :default)
@test isapprox(0.0, mean(X0_)[1]; atol=0.3)


# add the range measurement
addVariable!(fg, :x1, ContinuousEuclid{2})
eud = EuclidDistance(Normal(10,1))
fc = addFactor!(fg, [:x0;:x1], eud)
# X1_ = getBelief(fg, :x1) # not necessary really, can likely be deprecated

##

# project using utility function
X1_ = approxConvBelief(fg, fc, :x1)

##

@test isapprox(0, mean(X1_)[1], atol=1.5)
@test isapprox(0, mean(X1_)[2], atol=2.5) # looser for stochastic donut center

pts_ = getPoints(X1_)
@cast pts[i,j] := pts_[j][i]
N = size(pts, 2)

pts = collect(pts)
pts .^= 2
@test 0.5*N < sum( 7 .< sqrt.(sum(pts, dims=1)) .< 13 )

##

# current incremental solver builds a new tree and matches against old tree for recycling.
tree = IncrementalInference.buildTreeReset!(
  fg;
  drawpdf = false,
  show = false,
  ensureSolvable = false,
)

# setAllSolveFlags!(tree, false)

IncrementalInference.initTreeMessageChannels!(tree)
IncrementalInference.resetTreeCliquesForUpSolve!(tree)

##

retdict = IncrementalInference.upGibbsCliqueDensity(fg, tree.cliques[1], :default, LikelihoodMessage[])

## first tests with new HomotopyDensity works
# HomotopyDensity_legacy(retdict[:x0]) |> plotPoints
# HomotopyDensity_legacy(retdict[:x1]) |> plotPoints

##

cliq = 1
res = IncrementalInference.solveClique!(
  fg,
  tree,
  cliq
)

##
end


@testset "2D range-only donut posterior test, all up" begin

##

fg = initfg()

addVariable!(fg, :x0, ContinuousEuclid{2})
addFactor!(fg, [:x0], Prior(MvNormal(zeros(2),diagm([1;1.0]))))

addVariable!(fg, :x1, ContinuousEuclid{2})

eud = EuclidDistance(Normal(10,1))
addFactor!(fg, [:x0;:x1], eud)

##

tree = solveGraph!(fg)

##

@test isapprox(mean(getBelief(getState(fg, :x0, :default)))[1], 0, atol=1)
@test isapprox(mean(getBelief(getState(fg, :x0, :default)))[2], 0, atol=1)

pts_ = getBelief(fg, :x1) |> getPoints
@cast pts[i,j] := pts_[j][i]
N = size(pts, 2)

pts = collect(pts)
pts .^= 2
@test 0.5*N < sum( 7 .< sqrt.(sum(pts, dims=1)) .< 13 )

##
end


@testset "test upward clique message range density behavior" begin
## Test zero with on x and y-axis

N=100
points = [[100.0;0.0],[0.0;100.0]]
fg = IIF.generateGraph_EuclidDistance(points)

eo = [:x2; :x1; :l1]

fg_ = deepcopy(fg)
tree = buildTreeReset!(fg_, eo)

##

hist,upMessage = solveCliqUp!(fg_, tree, :x2; recordcliq=true);

##

sfg = hist[end].csmc.cliqSubFg
L1__ = getBelief(sfg, :l1) |> getPoints
@cast L1_[i,j] := L1__[j][i]

# check for for ring density

@test 0.2*N < sum( 0 .< L1_[1,:] .< 130)
@test 0.2*N < sum( -130 .< L1_[1,:] .< 0)
@test 0.2*N < sum( 100 .< L1_[2,:] .< 230)
@test 0.2*N < sum( -230 .< L1_[2,:] .< 100)

# and must be in a ring
L1_ = collect(L1_)
L1_[2,:] .-= 100
@test_broken 0.95*N < sum( 90 .< sqrt.(sum(L1_.^2, dims=1)) .< 110)

##

N=100
points = [[100.0;0.0],[0.0;100.0]]
fg = IIF.generateGraph_EuclidDistance(points)

# initVariable!(fg, :l1, [1000.0.*randn(2) for _ in 1:100])

## check regular full solution produces two modes


# similar test in RoME
for i in 1:1
  # global TP, N
  tree = solveGraph!(fg, eliminationOrder=eo);

  L1_ = getBelief(fg, :l1) |> getPoints
  @cast L1[i,j] := L1_[j][i] 
  # check that two modes exist
  @test (0.03*N < sum(-50 .< L1[1,:] .< 50))
  @test (0.03*N < sum(-50 .< L1[2,:] .< 50))
  # @error "suppressing dual mode tests, MUST restore before IIF v0.25, see #1305"
  @test (0.03*N < sum(50 .< L1[1,:] .< 150)) # always this one
  @test (0.03*N < sum(50 .< L1[2,:] .< 150))

end

# at least one of the 3 solves should produce the right result


##
#test one clique as in RoME
N=100
points = [[100.0;0.0],[0.0;100.0]]
fg = IIF.generateGraph_EuclidDistance(points)
getSolverParams(fg).graphinit = false

M = getManifold(fg, :l1)
TP = false
for i in 1:3
  # global TP, N
  tree = solveGraph!(fg);

  L1 = getBelief(fg, :l1) |> getPoints

  # check that two modes exist
  am1 = sum(isapprox.(Ref(M), L1, Ref([0.0,0.0]), atol=10))
  am2 = sum(isapprox.(Ref(M), L1, Ref([100.0,100.0]), atol=10))

  TP  = am1 > N*0.03
  TP &= am2 > N*0.03
  if TP 
    @info "test passed in $i"
    break
  end
end
@test_broken TP
##

end




@testset "Euclid Distance Tests" begin
# using Random
# Random.seed!(84)
# N=100
##
points = [[100.0],]
fg = IIF.generateGraph_EuclidDistance(points)
solveGraph!(fg)

@test isapprox(mean(getBelief(fg, :x1, :default))[1], 100, atol=1)

pts_ = getBelief(fg, :l1) |> getPoints
@cast pts[i,j] := pts_[j][i]
N = size(pts, 2)

# TODO add similar tests to the rest
@test_broken 0.3*N < sum(isapprox.(pts,  0, atol=5)) < 0.7*N
@test_broken 0.3*N < sum(isapprox.(pts,200, atol=5)) < 0.7*N

# Do it manually with a big inflation
# IIF._getCCW(fg, :x1l1f1).inflation = 100.0 # never gets there
# IIF._getCCW(fg, :x1l1f1).inflation = 150.0 # few iters gets there
fct = getObservation(fg, :x1l1f1)
deleteFactor!(fg, :x1l1f1)
addFactor!(fg, [:x1;:l1], fct; inflation=200.0)
# IIF._getCCW(fg, :x1l1f1).inflation = 200.0 # One almost, second good
pts = approxConv(fg, :x1l1f1, :l1)
initVariable!(fg, :l1, pts)
# plotKDE(fg, ls(fg))

pts_ = approxConv(fg, :x1l1f1, :l1)
initVariable!(fg, :l1, pts_)
# plotKDE(fg, ls(fg))

@cast pts[i,j] := pts_[j][i]

@test 0.3*N < sum(isapprox.(pts,  0, atol=5)) < 0.7*N
@test 0.3*N < sum(isapprox.(pts,200, atol=5)) < 0.7*N


## Test zero with x-axis
points = [[100.0;0.0],]
fg = IIF.generateGraph_EuclidDistance(points)
solveGraph!(fg)

## Test zero with y-axis
points = [[0.0;100.0],]
fg = IIF.generateGraph_EuclidDistance(points)
solveGraph!(fg)

## Test zero with xy-axis 2 points
points = [[0.0;100.0],[100.0;0.0]]
fg = IIF.generateGraph_EuclidDistance(points)
solveGraph!(fg)

## Test offsett with xy-axis 2 points
points = [[50.0;100.0],[100.0;50.0]]
fg = IIF.generateGraph_EuclidDistance(points; dist=50.0)
solveGraph!(fg)
# plotKDE(fg, ls(fg))

## Manual init
points = [[0.0;100.0],[100.0;0.0]]
fg = IIF.generateGraph_EuclidDistance(points)
getSolverParams(fg).inflation=3.0

initVariable!(fg, :x1, [rand(MvNormal([100.,0], [1.,1])) for _ in 1:N])
initVariable!(fg, :x2, [rand(MvNormal([0.,100], [1.,1])) for _ in 1:N])

# init = MixtureModel([MvNormal([100.,100], [10.,10]),
#                        MvNormal([0.,0], [10.,10])],
#                        [0.5, 0.5])
init = MvNormal([25.,25], [1.,1])
initVariable!(fg, :l1, [rand(init) for _ in 1:N])

# plotKDE(fg, ls(fg))

# normal 2 clique eliminationOrder 
eliminationOrder = [:l1; :x2; :x1]
# one clique eliminationOrder
eliminationOrder = [:l1; :x2; :x1]
tree = solveGraph!(fg; eliminationOrder)

##

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


##

# plotKDE(fg, :l1)

##

# pts = approxConv(fg, :x2l1f1, :l1)
# plotKDE(manikde!(ContinuousEuclid{2}, pts))
# plotLocalProduct(fg, :l1, levels=3)


## what would clique solution produce as up message


# @error "continue test dev with #1168"
#solve the clique in isolation
# hist = solveCliqUp!(fg, tree, :x2; recordcliq=true);
# printCliqHistorySummary(hist)
# sfg = hist[end].csmc.cliqSubFg

# the belief that would have been sent by this clique:
# L1 = IIF.getMessageBuffer(hist[11].csmc.cliq).upTx.belief[:l1] |> manikde!

# fnc_, csmc_ = repeatCSMStep!(hist, 5);
# sfg = csmc_.cliqSubFg

##

# plotKDE(sfg, :l1)

##

# initVariable!(sfg, :l1, pts)
# pts = approxConv(sfg, :x2l1f1, :l1)
# plotKDE(manikde!(ContinuousEuclid{2}, pts))