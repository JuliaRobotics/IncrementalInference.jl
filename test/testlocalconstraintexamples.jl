# test Odo, Obsv2

using IncrementalInference, KernelDensityEstimate
using Base.Test



fg = emptyFactorGraph()

N=100

doors = [0.0]'
pd = kde!(doors,[3.0])
pd = resample(pd,N);
bws = getBW(pd)[:,1]
doors2 = getPoints(pd);
v1 = addNode!(fg,:x1,doors,N=N)
f1  = addFactor!(fg,[v1],Obsv2( doors2, bws', [1.0])) #, samplefnc=getSample

tem = 2.0*randn(1,N)+getVal(v1)+50.0
v2 = addNode!(fg,:x2, tem, N=N)
odoc = Odo([50.0]',[2.0]',[1.0])
f2 = addFactor!(fg, [:x1; :x2], odoc ) #, samplefnc=getSample





# test evaluation of pose pose constraint
pts = evalFactor2(fg, f2, v2.index)
@test norm(Base.mean(pts,2)-[50.0]) < 10.0




# plotKDE( getVertKDE(fg,:x2) )














#
