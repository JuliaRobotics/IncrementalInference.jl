using IncrementalInference
# using KernelDensityEstimate, Gadfly # for vstack

fg = emptyFactorGraph()

N=100

doors = [-100.0;0.0;100.0;300.0]'
cov = [3.0]


v1 = addNode!(fg,:x1,ContinuousScalar,N=N)
f1  = addFactor!(fg,[v1], Obsv2(doors, cov', [1.0]))

tem = 2.0*randn(1,N)+getVal(v1)+50.0
v2 = addNode!(fg,:x2, ContinuousScalar, N=N)
addFactor!(fg,[v1;v2],Odo([50.0]',[2.0]',[1.0]))

v3=addNode!(fg,:x3,ContinuousScalar, N=N)
addFactor!(fg,[v2;v3],Odo([50.0]',[4.0]',[1.0]))

v4=addNode!(fg,:x4,ContinuousScalar, N=N)
addFactor!(fg,[v3;v4],Odo([50.0]',[2.0]',[1.0]))

v5=addNode!(fg,:x5,ContinuousScalar, N=N)
addFactor!(fg,[v4;v5],Odo([50.0]',[2.0]',[1.0]))


l1=addNode!(fg, :l1, ContinuousScalar, N=N)
addFactor!(fg, [v3,l1], Ranged([64.0],[0.5],[1.0]))

l2=addNode!(fg, :l2, ContinuousScalar, N=N)
addFactor!(fg, [v4,l2], Ranged([16.0],[0.5],[1.0]))

l3=addNode!(fg, :l3, ContinuousScalar, N=N)
addFactor!(fg, [v5,l3], Ranged([16.0],[0.5],[1.0]))


# writeGraphPdf(fg);

tree = prepBatchTree!(fg,drawpdf=true);


# do belief propagation inference over tree once
# inferOverTreeR!(fg, tree, N=100)
