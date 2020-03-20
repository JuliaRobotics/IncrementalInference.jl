
using Test
using IncrementalInference


@testset "basic test for tree initialization functionality" begin

# small canonical factor graph, without graphinit
fg = generateCanonicalFG_CaesarRing1D(graphinit=false)
getSolverParams(fg).graphinit = false

# do all init on tree as part of solve
tree, smt, hist = solveTree!(fg)


end


#
