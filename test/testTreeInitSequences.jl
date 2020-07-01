
using Test
using IncrementalInference


@testset "basic test for tree initialization functionality" begin

# small canonical factor graph, without graphinit
fg = generateCanonicalFG_CaesarRing1D(graphinit=false)
getSolverParams(fg).graphinit = false

# do all init on tree as part of solve
# getSolverParams(fg).drawtree = true
# tree, smt, hist = solveTree!(fg)


tree = wipeBuildNewTree!(fg)


getSolverParams(fg).limititers = 30

TS = Vector{Task}(undef,6)

TS[4] = solveCliq!(fg, tree, :x1, recordcliq=true, async=true)
TS[3] = solveCliq!(fg, tree, :l1, recordcliq=true, async=true)
TS[1] = solveCliq!(fg, tree, :x0, recordcliq=true, async=true)
TS[2] = solveCliq!(fg, tree, :x4, recordcliq=true, async=true)
TS[5] = solveCliq!(fg, tree, :x5, recordcliq=true, async=true)
TS[6] = solveCliq!(fg, tree, :x3, recordcliq=true, async=true)



end



using Logging





hist3

fnct3_7 = deepcopy(hist3[7][3])
csmc3_7 = deepcopy(hist3[7][4])

getCliqueData(csmc3_7.cliq).solveCondition = Condition()
csmc3_7.logger = SimpleLogger(stdout)


fnct3_8 = fnct3_7(csmc3_7)
csmc3_8 = deepcopy(csmc3_7)


Juno.@enter fnct3_8(csmc3_8)


#
