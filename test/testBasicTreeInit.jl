
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


# # debug only
# getSolverParams(fg).drawtree = true
# getSolverParams(fg).showtree = true
# getSolverParams(fg).dbg = true

# mkpath(getLogPath(fg))
# verbosefid = open(joinLogPath(fg, "csmVerbose.log"),"w")

tree, smt, hist = solveTree!(fg, timeout=100 ) # , verbose=true, verbosefid=verbosefid, timeout=50)

# flush(verbosefid)
# close(verbosefid)
# open(joinLogPath(fg, "csmLogicalReconstructMax.log"),"w") do io
#   IIF.reconstructCSMHistoryLogical(getLogPath(fg), fid=io)
# end


end


#
