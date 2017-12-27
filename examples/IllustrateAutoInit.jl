# example to illustate how autoinit functions work, and used as a development script towards a standard unit test.

using IncrementalInference


fg = emptyFactorGraph()

N=100

doors = reshape(Float64[-100.0;0.0;100.0;300.0],1,4)
pd = kde!(doors,[3.0])
pd = resample(pd,N);
bws = getBW(pd)[:,1]
doors2 = getPoints(pd);

v1 = addNode!(fg, :x0, autoinit=true)
f1  = addFactor!(fg,[v1],Obsv2( doors2, reshape(bws,1,1), [1.0])) #, samplefnc=getSample




#
