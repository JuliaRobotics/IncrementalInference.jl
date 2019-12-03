
using IncrementalInference
using Test
@testset "Test priors" begin

    N=200
    graphinits = [false, true]
    for  graphinit = graphinits

        fg = initfg()
        fg.solverParams.N = N


        addVariable!(fg, :x0, ContinuousScalar, autoinit = graphinit, N=N)
        addFactor!(fg, [:x0], Prior(Normal(-1.0, 1.0)))

        addVariable!(fg, :x1, ContinuousScalar, autoinit = graphinit, N=N)

        addVariable!(fg, :x2, ContinuousScalar, autoinit = graphinit, N=N)
        addFactor!(fg, [:x2], Prior(Normal(+1.0, 1.0)))

        addFactor!(fg, [:x0; :x1], LinearConditional(Normal(0.0, 0.1)))
        addFactor!(fg, [:x1; :x2], LinearConditional(Normal(0.0, 0.1)))

        #solve
        tree2, smt, hist = solveTree!(fg)
        x0_m = getKDEMean(getKDE(getVariable(fg, :x0)))[1]
        x1_m = getKDEMean(getKDE(getVariable(fg, :x1)))[1]
        x2_m = getKDEMean(getKDE(getVariable(fg, :x2)))[1]

        @info ("Testing means = 0 with 2 priors:\ngraphinit=$graphinit\nMeans: x0: $(x0_m), x1: $x1_m, x2: $x2_m")

        @test isapprox(x0_m, 0.0, atol = 0.1)
        @test isapprox(x1_m, 0.0, atol = 0.1)
        @test isapprox(x2_m, 0.0, atol = 0.1)

    end
end
