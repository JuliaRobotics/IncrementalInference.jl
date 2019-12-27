# test saving and loading of trees

using Test
using IncrementalInference

@testset "Test loading and saving of Bayes (Junction) tree" begin

    fg = loadCanonicalFG_Kaess(graphinit=false)
    tree = wipeBuildNewTree!(fg)

    # save and load tree
    saveTree(tree)
    tree2 = loadTree()

    # perform a few spot checks to see that the trees are similar
    @test length(tree.cliques) == length(tree2.cliques)
    @test getVariableOrder(tree) == getVariableOrder(tree2)

    for (clid,cl) in tree.cliques
      fsyms = getFrontals(cl)
      cl2 = getCliq(tree2, fsyms[1])
      fsyms2 = getFrontals(cl2)
      @test fsyms == fsyms2
      @test getCliqSeparatorVarIds(cl) == getCliqSeparatorVarIds(cl2)
      @test typeof(cl) == typeof(cl2)
    end

end
