using IncrementalInference


global fg = initfg()

global N=100

global doors = [-100.0;0.0;100.0;300.0]'
global cov = [3.0]


global v1 = addVariable!(fg,:x1,ContinuousScalar,N=N)
global f1  = addFactor!(fg,[v1], Prior(Normal(doors, cov)))

# tem = 2.0*randn(1,N)+getVal(v1)+50.0
global v2 = addVariable!(fg,:x2, ContinuousScalar, N=N)
addFactor!(fg,[v1;v2],Odo([50.0]',[2.0]',[1.0]))

global v3=addVariable!(fg,:x3,ContinuousScalar, N=N)
addFactor!(fg,[v2;v3],Odo([50.0]',[4.0]',[1.0]))

global v4=addVariable!(fg,:x4,ContinuousScalar, N=N)
addFactor!(fg,[v3;v4],Odo([50.0]',[2.0]',[1.0]))

global v5=addVariable!(fg,:x5,ContinuousScalar, N=N)
addFactor!(fg,[v4;v5],Odo([50.0]',[2.0]',[1.0]))


global l1=addVariable!(fg, :l1, ContinuousScalar, N=N)
addFactor!(fg, [v3,l1], Ranged([64.0],[0.5],[1.0]))

global l2=addVariable!(fg, :l2, ContinuousScalar, N=N)
addFactor!(fg, [v4,l2], Ranged([16.0],[0.5],[1.0]))

global l3=addVariable!(fg, :l3, ContinuousScalar, N=N)
addFactor!(fg, [v5,l3], Ranged([16.0],[0.5],[1.0]))


# writeGraphPdf(fg);

global tree = wipeBuildNewTree!(fg, drawpdf=true);
# global tree = prepBatchTree!(fg,drawpdf=true);

# run(`evince bt.pdf`)

# do belief propagation inference over tree once
# inferOverTreeR!(fg, tree, N=100)
