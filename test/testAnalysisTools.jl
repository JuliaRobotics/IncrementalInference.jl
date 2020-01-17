# Test for tree-based analysis tools found in `AnalysisTools.jl`.
using Test
using IncrementalInference

@testset "Number of non-zero calculation for frontals." begin
    # Alternative way of calculating number of upper triangular matrix elements.
    nnzFrontalsRecursive(dim) = dim==1 ? 1 : dim + nnzFrontalsRecursive(dim-1)
    # Both must agree for any integer number dimension.
    for dim in 1:100
        @test nnzFrontalsRecursive(dim) == nnzFrontals(dim)
    end
end

@testset "Number of non-zero calculation for full cliques." begin
    fg = generateCanonicalFG_Kaess()
    vo = [:l1, :l2, :x1, :x2, :x3]
    tree = buildTreeFromOrdering!(fg, vo)
    # Must agree with hand-calculated values, iSAM2 paper.
    @test nnzClique(tree.cliques[1]) == 3
    @test nnzClique(tree.cliques[2]) == 5
    @test nnzClique(tree.cliques[3]) == 2
end

@testset "Number of non-zero calculation for full trees." begin
    fg = generateCanonicalFG_Kaess()
    vo = [:l1, :l2, :x1, :x2, :x3]
    tree = buildTreeFromOrdering!(fg, vo)
    # Must agree with hand-calculated values, iSAM2 paper.
    @test nnzTree(tree) == 10
end
