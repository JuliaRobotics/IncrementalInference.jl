# test deconvolution functions

using Test
using IncrementalInference


@testset "Testing default deconvolution tools" begin

fg = generateCanonicalFG_CaesarRing1D()

tree, smt, hists = solveTree!(fg) 

msg = getMsgUpThis(tree.cliques[2])

tfg = addLikelihoodsDifferential!(msg)

# drawGraph(tfg, show=true)

@test intersect(ls(tfg), [:x2;:x6]) |> length == 2
@test lsf(tfg) |> length == 1
@test lsf(tfg, tags=[:UPWARD_DIFFERENTIAL]) |> length == 1

end





#
