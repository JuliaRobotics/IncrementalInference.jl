# test deconvolution functions

using Test
using IncrementalInference


##

@testset "basic deconvolution test" begin

##

fg = generateCanonicalFG_lineStep(2)

# drawGraph(fg, show=true)

## test trivial Prior

pred, meas = approxDeconv(fg, :x0f1)

@test mmd(pred, meas) < 1e-8

##

doautoinit!.(fg, [:x0; :x2])

##

pred, meas = approxDeconv(fg, :x0x2f1)

@test mmd(pred, meas) < 1e-5

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

X0 = getBelief(fg, :x0) |> getPoints

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

fg = generateCanonicalFG_CaesarRing1D()

# # TEMPORARY MUST COMMENT ON TRAVIS
# getSolverParams(fg).drawtree = true

mkpath(getLogPath(fg))
tree, smt, hists = solveTree!(fg, timeout=40, verbose=true) 

# msg = getMsgUpThis(tree.cliques[2])
msg = IIF.getMessageBuffer(tree.cliques[2]).upRx

tfg = buildCliqSubgraph(fg, tree.cliques[2])
addLikelihoodsDifferential!.(tfg, values(msg))

# drawGraph(tfg, show=true)

@test issetequal(ls(tfg), [:x2,:x4,:x6])
@test lsf(tfg) |> length == 2
@test lsf(tfg, tags=[:__UPWARD_DIFFERENTIAL__]) |> length == 2

##

end


@testset "deconv on <:AbstractRelativeMinimize" begin

##

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0], Prior(Normal()))
addFactor!(fg, [:x0;:x1], EuclidDistance(Normal(10,1)))

##

pts = approxConv(fg, :x0x1f1, :x1)

##

solveTree!(fg);

## make sure result is in the right place

@test abs(getPPE(fg, :x0).suggested[1]) < 1.0

X1 = getBelief(fg, :x1) |> getPoints
N = size(X1,2)
@test sum(-5 .< X1 .< 5) < 0.1*N
@test sum(X1 .< -15) < 0.1*N
@test sum(15 .< X1) < 0.1*N
@test 0.2*N .< sum(-15 .< X1 .< -5)
@test 0.2*N .< sum(5 .< X1 .< 15)

## not check deconv

pred, meas = approxDeconv(fg, :x0x1f1)

@test mmd(pred, meas) < 1e-2


##

# using KernelDensityEstimatePlotting, Gadfly
# Gadfly.set_default_plot_size(25cm,20cm)

# plotKDE([kde!(pred); kde!(meas)], c=["red"; "green"])


##

end


#
