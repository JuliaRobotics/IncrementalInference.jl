##


using IncrementalInference
using Test
using TensorCast
# using Statistics
# using Manifolds # should be done within regular exports

# using Pkg
# Pkg.resolve()

##

@testset "test Mixture sampling" begin

##

fg = initfg()
addVariable!(fg, :x0, ContinuousScalar)

mp = Mixture(Prior, (Normal(), Normal(10,1)),(1/2,1/2) )
addFactor!(fg, [:x0], mp)

##

pts_ = approxConv(fg, :x0f1, :x0)
@cast pts[i,j] := pts_[j][i]

N = size(pts,2)
@test 0.2*N < sum( -5 .< pts .< 5 )
@test 0.2*N < sum( 5 .< pts .< 15 )
@test sum( 15 .< pts  ) < 0.1*N
@test sum( pts .< -5 ) < 0.1*N
@test sum( 3 .< pts .< 7 ) < 0.1*N


# using KernelDensityEstimatePlotting, Gadfly
# Gadfly.set_default_plot_size(25cm,20cm)

# plotKDE(kde!(pts))

##

fg = initfg()
addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)

addFactor!(fg, [:x0], Prior(Normal()), graphinit=false)
initVariable!(fg, :x0, [zeros(1) for _ in 1:100])

mlr = Mixture(LinearRelative, (Normal(), Normal(10,1)),(1/2,1/2) )
addFactor!(fg, [:x0;:x1], mlr, graphinit=false)

##

pts_ = approxConv(fg, :x0x1f1, :x1)
@cast pts[i,j] := pts_[j][i]

# plotKDE(kde!(pts))

##

N = size(pts,2)
@test 0.2*N < sum( -5 .< pts .< 5 )
@test 0.2*N < sum( 5 .< pts .< 15 )
@test sum( 15 .< pts  ) < 0.1*N
@test sum( pts .< -5 ) < 0.1*N
@test sum( 3 .< pts .< 7 ) < 0.1*N


##

end



@testset "test packing of Mixture" begin

##

fg = initfg()
addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)

mp = Mixture(Prior, [Normal(); Normal(10,1)], [0.5;0.5])
f0 = addFactor!(fg, [:x0;], mp)

mr = Mixture(LinearRelative, (fancy=manikde!(ContinuousEuclid(1), [randn(1) for _ in 1:75]), naive=Normal(0,10)), [0.4;0.6])
f1 = addFactor!(fg, [:x0;:x1], mr)

##

pf0 = DFG.packFactor(f0)
pf1 = DFG.packFactor(f1)

# now test unpacking
fg_ = initfg();
addVariable!(fg_, :x0, ContinuousScalar)
addVariable!(fg_, :x1, ContinuousScalar)

##

f0_ = DFG.unpackFactor(fg_, pf0)
f1_ = DFG.unpackFactor(fg_, pf1)

##

# ENV["JULIA_DEBUG"] = "DistributedFactorGraphs"
@warn("Skipping pack/unpack compareFactor test for `timezone` and `zone`")
@show typeof(f1)
@show typeof(f1_)

@show  typeof(getSolverData(f1).fnc.varValsAll[]);
@show typeof(getSolverData(f1_).fnc.varValsAll[]);

@test DFG.compareFactor(f1, f1_, skip=[:components;:labels;:timezone;:zone;:vartypes;:fullvariables;:particleidx;:varidx])

@test IIF._getCCW(f1).usrfnc!.components.naive == IIF._getCCW(f1).usrfnc!.components.naive

# already ManifoldKernelDensity
A = IIF._getCCW(f1).usrfnc!.components.fancy
B = IIF._getCCW(f1_).usrfnc!.components.fancy

# A = ManifoldBelief(Euclid, IIF._getCCW(f1).usrfnc!.components.fancy )
# B = ManifoldBelief(Euclid, IIF._getCCW(f1_).usrfnc!.components.fancy )

@test mmd(A,B) < 1e-6

##

end


@testset "test simple Mixture" begin

##

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0], Prior(Normal(0.0,0.1)))

##

# require default ::UnitScaling constructor from all factors
mlr = Mixture(LinearRelative(I), [Normal(-1.0, 0.1), Normal(1.0, 0.1)], Categorical([0.5; 0.5]))

# test serialization while we are here
pmlr = convert(PackedMixture, mlr)
umlr = convert(Mixture, pmlr)

@test mlr.mechanics == umlr.mechanics
@test mlr.components == umlr.components
@test mlr.diversity == umlr.diversity

mlr = Mixture(LinearRelative, [Normal(-1.0, 0.1), Normal(1.0, 0.1)], Categorical([0.5; 0.5]))

addFactor!(fg, [:x0,:x1], mlr)

# To look at your factor graph
# if false
# using GraphPlot
# using DistributedFactorGraphs
# plotDFG(fg)
# end

tree = solveTree!(fg)

##

btd = getBelief(getVariable(fg, :x0))
@test isapprox(mean(getKDEfit(btd,distribution=Normal)), 0.0; atol=0.1) 

@test isapprox(std(getKDEfit(btd,distribution=Normal)), 0.1; atol=0.05) 

btd = getBelief(getVariable(fg, :x1))
pts_ = getPoints(btd)
@cast pts[i,j] := pts_[j][i]
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

##

end
