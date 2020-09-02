
using Test
using IncrementalInference

@testset "basic per clique stopping criteria" begin

fg = generateCanonicalFG_lineStep(1)
tree, smt, hist = solveTree!(fg, recordcliqs=[:x0;], limititercliqs=[(:x0=>2);])

@test haskey(hist, 1)

@test hist[1] |> length == 2

end


@testset "basic test for tree initialization functionality" begin

# small canonical factor graph, without graphinit
fg = generateCanonicalFG_CaesarRing1D(graphinit=false)
getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true
@show getLogPath(fg)

# temporary
# getSolverParams(fg).drawtree = true
# getSolverParams(fg).dbg = true

# do all init on tree as part of solve
# getSolverParams(fg).drawtree = true
tree, smt, hist = solveTree!(fg, verbose=true)


end


#
