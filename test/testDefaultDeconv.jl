# test deconvolution functions

using Test
using IncrementalInference
using TensorCast
using Manifolds: Euclidean

##

@testset "basic deconvolution test" begin

##

fg = generateGraph_LineStep(2)

# drawGraph(fg, show=true)

## test trivial Prior

pred, meas = approxDeconv(fg, :x0f1)

@test mmd(Euclidean(1),pred, meas) < 1e-8

##

doautoinit!.(fg, [:x0; :x2])

##

pred, meas = approxDeconv(fg, :x0x2f1)

@test mmd(Euclidean(1), pred, meas) < 1e-3

##

P_ = approxDeconvBelief(fg, :x0x2f1, LinearRelative)

@test isapprox( mean(Euclidean(1), meas), mean(P_), atol=0.2 )

##

end


# voodoo-lite
@testset "deconv through multihypo" begin

##

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :hypoA, ContinuousScalar)
addFactor!(fg, [:hypoA;], Prior(Normal(5,0.1)))
addVariable!(fg, :hypoB, ContinuousScalar)
addFactor!(fg, [:hypoB;], Prior(Normal(10,0.1)))

addFactor!(fg, [:x0; :hypoA; :hypoB], LinearRelative(Normal(10,0.1)), multihypo=[1;1/2;1/2])

solveTree!(fg);

##

# make sure each variable is where it should be first
@test isapprox(getPPE(fg, :hypoA).suggested[1], 5, atol=1)
@test isapprox(getPPE(fg, :hypoB).suggested[1], 10,atol=1)

X0_ = getBelief(fg, :x0)
X0 = AMP._pointsToMatrixCoords(X0_.manifold, getPoints(X0_))
# TensorCast.@cast X0[i,j] := X0_[j][i]

N = size(X0,2)
@test 0.2*N < sum( -7.5 .< X0 .< -2.5 )
@test 0.2*N < sum( -2.5 .< X0 .< 2.5 )
@test sum( 2.5 .< X0 ) < 0.05*N
@test sum( X0 .< -7.5 ) < 0.05*N
@test sum( -3.5 .< X0 .< -1.5 ) < 0.1*N

## do deconv and check

@error "approxDeconv on multihypo not fixed yet, see #467, #927"
# pred, meas = approxDeconv(fg, lsf(fg, LinearRelative)[1])

##

end


@testset "deconvolution tools via differential factors" begin

##

fg = generateGraph_CaesarRing1D()
getSolverParams(fg).useMsgLikelihoods = true

vo = [:x3,:x5,:x1,:l1,:x4,:x2,:x6,:x0]

mkpath(getLogPath(fg))
tree = solveTree!(fg, eliminationOrder=vo, verbose=true) #, timeout=5) # timeout creates interrupt exception

msg = IIF.getMessageBuffer(getClique(tree,2)).upRx

tfg = buildCliqSubgraph(fg, getClique(tree,2))
addLikelihoodsDifferential!.(tfg, values(msg))

# drawGraph(tfg, show=true)

@show ls(tfg)
@test issetequal(ls(tfg), [:x2,:x4,:x6])
@test lsf(tfg) |> length == 2
@test lsf(tfg, tags=[:__UPWARD_DIFFERENTIAL__]) |> length == 2

##

end


@testset "deconv on <:AbstractRelativeMinimize" begin

##

fg = initfg()
getSolverParams(fg).useMsgLikelihoods = true

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0], Prior(Normal()))
doautoinit!(fg,:x0)
addFactor!(fg, [:x0;:x1], EuclidDistance(Normal(10,1)))

##

# initAll!(fg)

pts = approxConv(fg, :x0x1f1, :x1)

##

solveTree!(fg);

## make sure result is in the right place

@test abs(getPPE(fg, :x0).suggested[1]) < 1.0

X1_ = getBelief(fg, :x1) |> getPoints
TensorCast.@cast X1[i,j] := X1_[j][i]

N = size(X1,2)
@test sum(-5 .< X1 .< 5) < 0.1*N
@test sum(X1 .< -15) < 0.1*N
@test sum(15 .< X1) < 0.1*N
@test 0.2*N .< sum(-15 .< X1 .< -5)
@test 0.2*N .< sum(5 .< X1 .< 15)

## not check deconv

pred, meas = approxDeconv(fg, :x0x1f1)

@test mmd(Euclidean(1), pred, meas) < 1e-1


##

# using KernelDensityEstimatePlotting, Gadfly
# Gadfly.set_default_plot_size(25cm,20cm)

# plotKDE([kde!(pred); kde!(meas)], c=["red"; "green"])


##

end



#
