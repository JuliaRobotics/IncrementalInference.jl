
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
getSolverParams(fg).N = N
getSolverParams(fg).graphinit = graphinit
getSolverParams(fg).treeinit = !graphinit

  addVariable!(fg, :x0, ContinuousScalar, N=N)
  addFactor!(fg, [:x0], Prior(Normal(-1.0, 1.0)))

  addVariable!(fg, :x1, ContinuousScalar, N=N)

  addVariable!(fg, :x2, ContinuousScalar, N=N)
  addFactor!(fg, [:x2], Prior(Normal(+1.0, 1.0)))

  addFactor!(fg, [:x0; :x1], LinearRelative(Normal(0.0, 0.01)))
  addFactor!(fg, [:x1; :x2], LinearRelative(Normal(0.0, 0.01)))

  #solve
  tree = solveTree!(fg)
  x0_m = mean(getBelief(getState(fg, :x0, :default)))
  x1_m = mean(getBelief(getState(fg, :x1, :default)))
  x2_m = mean(getBelief(getState(fg, :x2, :default)))

  @info ("Testing means = 0 with 2 priors:\ngraphinit=$graphinit\nMeans: x0: $(x0_m), x1: $x1_m, x2: $x2_m")

  @test_skip isapprox(x0_m, 0.0, atol = 0.1)
  @test_skip isapprox(x1_m, 0.0, atol = 0.1)
  @test_skip isapprox(x2_m, 0.0, atol = 0.1)

  @warn "priorusetest.jl is testing with large tolerances"
  @test isapprox(x0_m[1], 0.0, atol = 1.0)
  @test isapprox(x1_m[1], 0.0, atol = 1.25)
  @test isapprox(x2_m[1], 0.0, atol = 1.5)

  #testing if values are close to one another
  testvals = [x0_m, x1_m, x2_m]
  meanval = mean(testvals)
  @test_skip all(isapprox.(testvals[1], meanval, atol=0.4))

end

##

for  graphinit = graphinits

fg = initfg()
getSolverParams(fg).N = N
getSolverParams(fg).graphinit = graphinit
getSolverParams(fg).treeinit = !graphinit

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

  x0_m = mean(getBelief(getState(fg, :x0, :default)))
  x1_m = mean(getBelief(getState(fg, :x1, :default)))
  x2_m = mean(getBelief(getState(fg, :x2, :default)))
  l0_m = mean(getBelief(getState(fg, :l0, :default)))
  l1_m = mean(getBelief(getState(fg, :l1, :default)))

  @info ("Testing means = 0 with 2 priors:\ngraphinit=$graphinit\nMeans: x0: $(x0_m), x1: $x1_m, x2: $x2_m, l0: $l0_m, l1: $l1_m")

  @test_skip isapprox(x0_m[1], 0.0, atol = 0.1)
  @test_skip isapprox(x1_m[1], 0.0, atol = 0.1)
  @test_skip isapprox(x2_m[1], 0.0, atol = 0.1) #1.2)
  @test_skip isapprox(l0_m[1], 0.0, atol = 0.1) #1.5)
  @test_skip isapprox(l1_m[1], 0.0, atol = 0.1) #1.2)

  #testing if values are close to one another
  @show testvals = [x0_m, x1_m, x2_m, l0_m, l1_m]
  @show meanval = mean(testvals)
  @test isapprox(testvals[1], meanval, atol=0.3) # skip?
  @test isapprox(testvals[2], meanval, atol=0.3) # skip?
  @test isapprox(testvals[3], meanval, atol=0.3) # skip?
  @test isapprox(testvals[4], meanval, atol=0.3) # skip?
  @test isapprox(testvals[5], meanval, atol=0.3) # skip?

end

##

end
