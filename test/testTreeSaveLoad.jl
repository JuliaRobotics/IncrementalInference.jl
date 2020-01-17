# test saving and loading of trees

using Test
using IncrementalInference

@testset "Test loading and saving of Bayes (Junction) tree" begin

    fg = generateCanonicalFG_Kaess(graphinit=false)
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


@testset "Test loading and saving of Bayes (Junction) tree" begin

    fg = generateCanonicalFG_Kaess(graphinit=false)
    tree = wipeBuildNewTree!(fg)

    # save and load tree as array
    filepath = saveTree([tree;deepcopy(tree)])
    trees = loadTree(filepath)

    # perform a few spot checks to see that the trees are similar
    @test length(tree.cliques) == length(trees[1].cliques)
    @test getVariableOrder(tree) == getVariableOrder(trees[1])

    for (clid,cl) in tree.cliques
      fsyms = getFrontals(cl)
      cl2 = getCliq(trees[1], fsyms[1])
      fsyms2 = getFrontals(cl2)
      @test fsyms == fsyms2
      @test getCliqSeparatorVarIds(cl) == getCliqSeparatorVarIds(cl2)
      @test typeof(cl) == typeof(cl2)
    end

end
