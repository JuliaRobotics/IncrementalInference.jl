# test saving and loading of trees

using Test
using IncrementalInference
using MetaGraphs
using Graphs

##

@testset "Test loading and saving of Bayes (Junction) tree" begin
##

fg = generateGraph_Kaess(graphinit=false)
tree = buildTreeReset!(fg)

# save and load tree
saveTree(tree)
tree2 = loadTree()

# perform a few spot checks to see that the trees are similar
@test length(tree.cliques) == length(tree2.cliques)
@test getEliminationOrder(tree) == getEliminationOrder(tree2)

for (clid,cl) in tree.cliques
  fsyms = getFrontals(cl)
  cl2 = getClique(tree2, fsyms[1])
  fsyms2 = getFrontals(cl2)
  @test fsyms == fsyms2
  @test getCliqSeparatorVarIds(cl) == getCliqSeparatorVarIds(cl2)
  @test typeof(cl) == typeof(cl2)
end

##
end


@testset "Test loading and saving of Bayes (Junction) tree" begin
##

fg = generateGraph_Kaess(graphinit=false)
tree = buildTreeReset!(fg)

# save and load tree as array
filepath = saveTree([tree;deepcopy(tree)])
trees = loadTree(filepath)

# perform a few spot checks to see that the trees are similar
@test length(tree.cliques) == length(trees[1].cliques)
@test getEliminationOrder(tree) == getEliminationOrder(trees[1])

for (clid,cl) in tree.cliques
  fsyms = getFrontals(cl)
  cl2 = getClique(trees[1], fsyms[1])
  fsyms2 = getFrontals(cl2)
  @test fsyms == fsyms2
  @test getCliqSeparatorVarIds(cl) == getCliqSeparatorVarIds(cl2)
  @test typeof(cl) == typeof(cl2)
end

##
end
