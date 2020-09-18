# test deconvolution functions

using Test
using IncrementalInference


@testset "Testing default deconvolution tools" begin

fg = generateCanonicalFG_CaesarRing1D()

# # TEMPORARY MUST COMMENT ON TRAVIS
# getSolverParams(fg).drawtree = true

mkpath(getLogPath(fg))
verbosefid = open(joinLogPath(fg, "csmVerbose.log"),"w")
tree, smt, hists = solveTree!(fg, timeout=40, verbose=true, verbosefid=verbosefid) 
flush(verbosefid)
close(verbosefid)
open(joinLogPath(fg, "csmLogicalReconstructMax.log"),"w") do io
  IIF.reconstructCSMHistoryLogical(getLogPath(fg), fid=io)
end

msg = getMsgUpThis(tree.cliques[2])

tfg = addLikelihoodsDifferential!(msg)

# drawGraph(tfg, show=true)

@test intersect(ls(tfg), [:x2;:x6]) |> length == 2
@test lsf(tfg) |> length == 1
@test lsf(tfg, tags=[:UPWARD_DIFFERENTIAL]) |> length == 1

end





#
