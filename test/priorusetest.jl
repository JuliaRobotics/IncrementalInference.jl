
# using Revise

using IncrementalInference
using Test

##

@testset "Test priors" begin

##

N=100
graphinits = [false, true]

# TEMP value when not using for loop
graphinit = true

for  graphinit = graphinits

fg = initfg()
fg.solverParams.N = N
fg.solverParams.graphinit = graphinit
fg.solverParams.treeinit = !graphinit

addVariable!(fg, :x0, ContinuousScalar, N=N)
addFactor!(fg, [:x0], Prior(Normal(-1.0, 1.0)))

addVariable!(fg, :x1, ContinuousScalar, N=N)

addVariable!(fg, :x2, ContinuousScalar, N=N)
addFactor!(fg, [:x2], Prior(Normal(+1.0, 1.0)))

addFactor!(fg, [:x0; :x1], LinearRelative(Normal(0.0, 0.01)))
addFactor!(fg, [:x1; :x2], LinearRelative(Normal(0.0, 0.01)))

#solve
tree = solveTree!(fg)
x0_m = getKDEMean(getBelief(getVariable(fg, :x0)))[1]
x1_m = getKDEMean(getBelief(getVariable(fg, :x1)))[1]
x2_m = getKDEMean(getBelief(getVariable(fg, :x2)))[1]

@info ("Testing means = 0 with 2 priors:\ngraphinit=$graphinit\nMeans: x0: $(x0_m), x1: $x1_m, x2: $x2_m")

@test_skip isapprox(x0_m, 0.0, atol = 0.1)
@test_skip isapprox(x1_m, 0.0, atol = 0.1)
@test_skip isapprox(x2_m, 0.0, atol = 0.1)

@warn "priorusetest.jl is testing with large tolerances"
@test isapprox(x0_m, 0.0, atol = 1.0)
@test isapprox(x1_m, 0.0, atol = 1.0)
@test isapprox(x2_m, 0.0, atol = 1.0)

#testing if values are close to one another
testvals = [x0_m, x1_m, x2_m]
meanval = mean(testvals)
@test all(isapprox.(testvals, meanval, atol=0.4))

end

##

for  graphinit = graphinits

fg = initfg()
fg.solverParams.N = N
fg.solverParams.graphinit = graphinit
fg.solverParams.treeinit = !graphinit

addVariable!(fg, :x0, ContinuousScalar, N=N)
addFactor!(fg, [:x0], Prior(Normal(-1.0, 1.0)))

addVariable!(fg, :l0, ContinuousScalar, N=N)
addFactor!(fg, [:l0], Prior(Normal(+1.0, 1.0)))

addVariable!(fg, :l1, ContinuousScalar, N=N)

addFactor!(fg, [:x0; :l0], LinearRelative(Normal(0, 0.01)))
addFactor!(fg, [:x0; :l1], LinearRelative(Normal(0, 0.01)))

addVariable!(fg, :x1, ContinuousScalar, N=N)
addFactor!(fg, [:x0; :x1], LinearRelative(Normal(0, 0.01)))

addVariable!(fg, :x2, ContinuousScalar, N=N)
addFactor!(fg, [:x1; :x2], LinearRelative(Normal(0, 0.01)))

addFactor!(fg, [:x2; :l0], LinearRelative(Normal(0, 0.01)))
addFactor!(fg, [:x2; :l1], LinearRelative(Normal(0, 0.01)))

#solve
tree = solveTree!(fg)

x0_m = getKDEMean(getBelief(getVariable(fg, :x0)))[1]
x1_m = getKDEMean(getBelief(getVariable(fg, :x1)))[1]
x2_m = getKDEMean(getBelief(getVariable(fg, :x2)))[1]
l0_m = getKDEMean(getBelief(getVariable(fg, :l0)))[1]
l1_m = getKDEMean(getBelief(getVariable(fg, :l1)))[1]

@info ("Testing means = 0 with 2 priors:\ngraphinit=$graphinit\nMeans: x0: $(x0_m), x1: $x1_m, x2: $x2_m, l0: $l0_m, l1: $l1_m")

@test_skip isapprox(x0_m, 0.0, atol = 0.1)
@test_skip isapprox(x1_m, 0.0, atol = 0.1)
@test_skip isapprox(x2_m, 0.0, atol = 0.1)
@test_skip isapprox(l0_m, 0.0, atol = 0.1)
@test_skip isapprox(l1_m, 0.0, atol = 0.1)

@warn "priorusetest.jl is testing with large tolerances"
@test isapprox(x0_m, 0.0, atol = 1.0)
@test isapprox(x1_m, 0.0, atol = 1.0)
@test isapprox(x2_m, 0.0, atol = 1.0)
@test isapprox(l0_m, 0.0, atol = 1.2)
@test isapprox(l1_m, 0.0, atol = 1.2)

#testing if values are close to one another
@show testvals = [x0_m, x1_m, x2_m, l0_m, l1_m]
@show meanval = mean(testvals)
@test all(isapprox.(testvals, meanval, atol=0.3))

end

##

end
