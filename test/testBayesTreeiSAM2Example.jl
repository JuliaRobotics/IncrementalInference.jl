using IncrementalInference
using Test

# import IncrementalInference: buildCliquePotentials


global N=100
global fg = initfg()

# doors = [-100.0;0.0;100.0;300.0]'
# cov = [3.0]

addVariable!(fg,:x1, ContinuousScalar, N=N)
addFactor!(fg, [:x1;], Prior(Normal()))


addVariable!(fg,:x2, ContinuousScalar, N=N)
addFactor!(fg,[:x1, :x2], LinearConditional(Normal()))

addVariable!(fg, :x3, ContinuousScalar, N=N)
addFactor!(fg,[:x2,:x3],LinearConditional(Normal()))


addVariable!(fg, :l1, ContinuousScalar, N=N)
addFactor!(fg, [:x1,:l1], LinearConditional(Normal()) )
addFactor!(fg, [:x2,:l1], LinearConditional(Normal()) )

addVariable!(fg, :l2, ContinuousScalar, N=N)
addFactor!(fg, [:x3,:l2], LinearConditional(Normal()))


@testset "test building tree native" begin

global fg

resetFactorGraphNewTree!(fg)

# p = getEliminationOrder(fg, ordering=:qr)
p = [:l1;:l2;:x1;:x2;:x3]

println()
global fge = deepcopy(fg)
println("Building Bayes net...")
buildBayesNet!(fge, p)

global tree = emptyBayesTree()
buildTree!(tree, fge, p)


@test num_vertices(tree.bt) == 3


# Michael reference -- x2->x1, x2->x3, x2->x4, x2->l1, x4->x3, l1->x3, l1->x4
# if false
#   println("Bayes Tree")
#   # Graphs.plot(tree.bt)
#   fid = open("bt.dot","w+")
#   write(fid,Graphs.to_dot(tree.bt))
#   close(fid)
#   run(`dot bt.dot -Tpdf -o bt.pdf`)
# end

println("Find potential functions for each clique")
cliq = tree.cliques[1] # start at the root
buildCliquePotentials(fg, tree, cliq); # fg does not have the marginals as fge does

end



# TODO -- add testing to ensure this is the correct tree!

@warn "add test tree verification"

# run(`evince bt.pdf`)


@testset "build tree from ordering" begin

global fg

resetFactorGraphNewTree!(fg)
vo = getEliminationOrder(fg)
tree = buildTreeFromOrdering!(fg, vo)

end



#
