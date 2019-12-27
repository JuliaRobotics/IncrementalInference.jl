# another tree test

using IncrementalInference
using Test


@testset "Variable ordering Bayes tree member check." begin

  fg = loadCanonicalFG_Kaess()
  # Choose specific variable ordering and perform check.
  vo = [:l1, :l2, :x1, :x2, :x3]
  tree = buildTreeFromOrdering!(fg, vo)
  @test vo == tree.variableOrder
  @test vo == getVariableOrder(tree)
  @test vo == getEliminationOrder(tree)

end


@testset "Test tree formation and consistent APIs" begin

  fg = loadCanonicalFG_TestSymbolic()

  #writeGraphPdf(fg, show=true)

  eo = [:x1; :l3; :l1; :x5; :x2; :l2; :x4; :x3]

  tree = buildTreeFromOrdering!(fg,eo)
  # drawTree(tree, show=true)

  @warn "TODO, complete test on tree formation"


  ## test variable order APIs consistent, see issue 499

  vo = getEliminationOrder(fg)
  tree1 = resetBuildTreeFromOrder!(fg, vo)
  # drawTree(tree1, show=true)

  tree2 = wipeBuildNewTree!(fg)
  # drawTree(tree2, show=true)

  @test getVariableOrder(tree1) == getVariableOrder(tree1)

end






#
