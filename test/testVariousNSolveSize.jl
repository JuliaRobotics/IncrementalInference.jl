# test for variations in DOF solving,  #227 #316 #430

using IncrementalInference
using Test


@testset "Test for variations in N while solving" begin

fg = loadCanonicalFG_CaesarRing1D(graphinit=true)

# 150 is non-standard
getSolverParams(fg).N = 150
tree, smt, hist = solveTree!(fg)

# test with change for second solve
getSolverParams(fg).N = 200
tree, smt, hist = solveTree!(fg)

end

#
