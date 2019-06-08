using IncrementalInference
using Graphs
# using KernelDensityEstimate, Gadfly # for vstack

using Test


global N=100
global fg = initfg()

# doors = [-100.0;0.0;100.0;300.0]'
# cov = [3.0]

global v1 = addVariable!(fg,:x1, ContinuousScalar, N=N)
global f1  = addFactor!(fg, [:x1;], Prior(Normal()))

# tem = 2.0*randn(1,N)+getVal(v1)+50.0
global v2 = addVariable!(fg,:x2, ContinuousScalar, N=N)
addFactor!(fg,[:x1, :x2], LinearConditional(Normal()))

global v3=addVariable!(fg, :x3, ContinuousScalar, N=N) # 4.0*randn(1,N)+getVal(v2)+50.0
addFactor!(fg,[:x2,:x3],LinearConditional(Normal()))


global l1=addVariable!(fg, :l1, ContinuousScalar, N=N) # 0.5*randn(1,N)+getVal(v3)+64.0
addFactor!(fg, [:x1,:l1], LinearConditional(Normal()) )
addFactor!(fg, [:x2,:l1], LinearConditional(Normal()) )

global l2=addVariable!(fg, :l2, ContinuousScalar, N=N) # 0.5*randn(1,N)+getVal(v3)+64.0
addFactor!(fg, [:x3,:l2], LinearConditional(Normal()))

# addVariable!(fg, :x4, ContinuousScalar, N=N) # 4.0*randn(1,N)+getVal(v2)+50.0
# addFactor!(fg,[:x3,:x4],LinearConditional(Normal()))




# for thesis
# v4=addVariable!(fg,:x4,4.0*randn(1,N)+getVal(v2)+50.0, N=N)
# addFactor!(fg,[:x4,:l2],Odo([50.0]',[4.0]',[1.0]))



# Graphs.plot(fg.g)
# writeGraphPdf(fg);
# run(`evince fg.pdf`)


# p = IncrementalInference.getEliminationOrder(fg, ordering=:qr)
# global p = [7,10,1,3,5,12];
# global p = [7,10,1,3,5];
p = [:l1;:l2;:x1;:x2;:x3]

println()
global fge = deepcopy(fg)
println("Building Bayes net...")
buildBayesNet!(fge, p)

global tree = emptyBayesTree()
buildTree!(tree, fge, p)

# println("Bayes Net")
# sleep(0.1)
#fid = open("bn.dot","w+")
#write(fid,to_dot(fge.bn))
#close(fid)


@test num_vertices(tree.bt) == 3


# Michael reference -- x2->x1, x2->x3, x2->x4, x2->l1, x4->x3, l1->x3, l1->x4
if false
  println("Bayes Tree")
  # Graphs.plot(tree.bt)
  fid = open("bt.dot","w+")
  write(fid,Graphs.to_dot(tree.bt))
  close(fid)
  run(`dot bt.dot -Tpdf -o bt.pdf`)
end

# drawTree(tree)

import IncrementalInference: buildCliquePotentials

println("Find potential functions for each clique")
cliq = tree.cliques[1] # start at the root
buildCliquePotentials(fg, tree, cliq); # fg does not have the marginals as fge does




# TODO -- add testing to ensure this is the correct tree!

@warn "add test tree verification"

# run(`evince bt.pdf`)





#
