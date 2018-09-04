using IncrementalInference
using Graphs
# using KernelDensityEstimate, Gadfly # for vstack

using Base: Test

fg = emptyFactorGraph()

N=100

doors = [-100.0;0.0;100.0;300.0]'
cov = [3.0]


v1 = addNode!(fg,:x1, ContinuousScalar, N=N)
f1  = addFactor!(fg, [:x1;], Obsv2(doors, cov', [1.0]))

tem = 2.0*randn(1,N)+getVal(v1)+50.0
v2 = addNode!(fg,:x2, ContinuousScalar, N=N)
addFactor!(fg,[:x1, :x2],Odo([50.0]',[2.0]',[1.0]))

v3=addNode!(fg, :x3, ContinuousScalar, N=N) # 4.0*randn(1,N)+getVal(v2)+50.0
addFactor!(fg,[:x2,:x3],Odo([50.0]',[4.0]',[1.0]))


l1=addNode!(fg, :l1, ContinuousScalar, N=N) # 0.5*randn(1,N)+getVal(v3)+64.0
addFactor!(fg, [:x1,:l1], Ranged([64.0],[0.5],[1.0]))
addFactor!(fg, [:x2,:l1], Ranged([16.0],[0.5],[1.0]))

l2=addNode!(fg, :l2, ContinuousScalar, N=N) # 0.5*randn(1,N)+getVal(v3)+64.0
addFactor!(fg, [:x3,:l2], Ranged([16.0],[0.5],[1.0]))



# for thesis
# v4=addNode!(fg,:x4,4.0*randn(1,N)+getVal(v2)+50.0, N=N)
# addFactor!(fg,[:x4,:l2],Odo([50.0]',[4.0]',[1.0]))



# Graphs.plot(fg.g)
# writeGraphPdf(fg);
# run(`evince fg.pdf`)


# p = IncrementalInference.getEliminationOrder(fg, ordering=:qr)
p = [7,10,1,3,5,12];
p = [7,10,1,3,5];

println()
fge = deepcopy(fg)
println("Building Bayes net...")
buildBayesNet!(fge, p)

tree = emptyBayesTree()
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

# GraphViz.Graph(to_dot(tree.bt))
#Michael reference 3sig -- x2l1x4x3    x1|x2

import IncrementalInference: buildCliquePotentials

println("Find potential functions for each clique")
cliq = tree.cliques[1] # start at the root
buildCliquePotentials(fg, tree, cliq); # fg does not have the marginals as fge does

# now update all factor graph vertices used for this tree
for v in Graphs.vertices(fg.g)
  IncrementalInference.dlapi.updatevertex!(fg, v)
end




# TODO -- add testing to ensure this is the correct tree!



# run(`evince bt.pdf`)
