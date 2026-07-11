# another tree test


using IncrementalInference
using Test


@testset "Variable ordering Bayes tree member check." begin

  fg = generateGraph_Kaess()
  # Choose specific variable ordering and perform check.
  vo = [:l1, :l2, :x1, :x2, :x3]
  tree = buildTreeReset!(fg, vo)
  @test vo == tree.eliminationOrder
  @test vo == getEliminationOrder(tree)

end

@testset "Test Caesar Ring 1D symbolic tree construction" begin

fg = generateGraph_CaesarRing1D()
# drawGraph(fg, show=true)

eo = [:x0;:x2;:x4;:x6;:x1;:l1;:x5;:x3;]

tree = buildTreeReset!(fg,eo)
# drawTree(tree, show=true)

@test length(tree.cliques) == 6

C0 = getClique(tree, :x3)
@test intersect( C0 |> getFrontals, [:x3; :x5; :l1]) |> length == 3
@test C0 |> getCliqSeparatorVarIds |> length == 0

cC0 = getChildren(tree, C0)
@test cC0 |> length == 3

C1 = getClique(tree, :x1)
@test C1 in cC0
@test C1 |> getFrontals == [:x1;]
@test intersect(C1 |> getCliqSeparatorVarIds, [:x3;:l1]) |> length == 2

cC1 = getChildren(tree, C1)
@test cC1 |> length == 2

C4 = getClique(tree, :x2)
@test C4 in cC1
@test C4 |> getFrontals == [:x2;]
@test intersect(C4 |> getCliqSeparatorVarIds, [:x3;:x1]) |> length == 2

C5 = getClique(tree, :x0)
@test C5 in cC1
@test C5 |> getFrontals == [:x0;]
@test intersect(C5 |> getCliqSeparatorVarIds, [:l1;:x1]) |> length == 2

C2 = getClique(tree, :x6)
@test C2 in cC0
@test C2 |> getFrontals == [:x6;]
@test intersect(C2 |> getCliqSeparatorVarIds, [:l1;:x5]) |> length == 2

C3 = getClique(tree, :x4)
@test C3 in cC0
@test C3 |> getFrontals == [:x4;]
@test intersect(C3 |> getCliqSeparatorVarIds, [:x3;:x5]) |> length == 2

end


@testset "Test tree formation and consistent APIs" begin

fg = generateGraph_TestSymbolic()

#writeGraphPdf(fg, show=true)

eo = [:x1; :l3; :l1; :x5; :x2; :l2; :x4; :x3]

tree = buildTreeReset!(fg,eo)
# drawTree(tree, show=true)

@warn "TODO, complete further testing on tree formation"


## test variable order APIs consistent, see issue 499

vo = getEliminationOrder(fg)
tree1 = buildTreeReset!(fg, vo)
# drawTree(tree1, show=true)

tree2 = buildTreeReset!(fg)
# drawTree(tree2, show=true)

@test getEliminationOrder(tree1) == getEliminationOrder(tree1)

end






#
