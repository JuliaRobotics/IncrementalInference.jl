using IncrementalInference, KernelDensityEstimate

fg = emptyFactorGraph()

N=100

doors = reshape(Float64[-100.0;0.0;100.0;300.0],1,4)
pd = kde!(doors,[3.0])
pd = resample(pd,N);
bws = getBW(pd)[:,1]
doors2 = getPoints(pd);

v1 = addNode!(fg,:x0,ContinuousScalar,N=N)
f1  = addFactor!(fg,[v1],Obsv2( doors2, reshape(bws,1,1), [1.0])) #, samplefnc=getSample

# tem = 2.0*randn(1,N)+getVal(v1)+50.0
v2 = addNode!(fg,:x2, ContinuousScalar, N=N)
addFactor!(fg, [:x0; :x2], Odo(50.0*ones(1,1),2.0*ones(1,1),[1.0])) #, samplefnc=getSample
# addFactor!(fg, [v1;v2], Odo(50.0*ones(1,1),[2.0]',[1.0])) #, samplefnc=getSample


# monocular sighting would look something like
#addFactor!(fg, Mono, [:x3,:l1], [14.0], [1.0], [1.0])
#addFactor!(fg, Mono, [:x4,:l1], [11.0], [1.0], [1.0])

v3=addNode!(fg,:x3,ContinuousScalar, N=N)
addFactor!(fg,[v2;v3],Odo(50.0*ones(1,1),4.0*ones(1,1),[1.0])) #, samplefnc=getSample
f2 = addFactor!(fg,[v3], Obsv2(doors2, bws', [1.0]))

v4=addNode!(fg,:x4,ContinuousScalar, N=N)
addFactor!(fg,[v3;v4],Odo(50.0*ones(1,1),2.0*ones(1,1),[1.0])) #, samplefnc=getSample


l1=addNode!(fg, :l1, ContinuousScalar, N=N)
addFactor!(fg, [v3,l1], Ranged([64.0],[0.5],[1.0])) #, samplefnc=getSample
addFactor!(fg, [v4,l1], Ranged([16.0],[0.5],[1.0])) #, samplefnc=getSample



v5=addNode!(fg,:x5,ContinuousScalar, N=N)
addFactor!(fg,[v4;v5],Odo(50.0*ones(1,1),2.0*ones(1,1),[1.0])) #, samplefnc=getSample


v6=addNode!(fg,:x6,ContinuousScalar, N=N)
addFactor!(fg,[v5;v6],Odo(40.0*ones(1,1),1.20*ones(1,1),[1.0])) #, samplefnc=getSample


v7=addNode!(fg,:x7,ContinuousScalar, N=N)
addFactor!(fg,[v6;v7],Odo(60.0*ones(1,1),2.0*ones(1,1),[1.0])) #, samplefnc=getSample

f3 = addFactor!(fg,[v7], Obsv2(doors, reshape(bws,1,1), [1.0])) #, samplefnc=getSample


# HMM computed ground truth, extended for 7 poses with landmark
gt = Dict{Symbol, Array{Float64,2}}()
gt[:x0]=reshape(Float64[0.0;1.97304 ],2,1) # -0.0342366
gt[:x2]=reshape(Float64[50.0; 2.83153 ],2,1) # 49.8797
gt[:x3]=reshape(Float64[100.0; 1.65557 ],2,1) # 99.8351
gt[:x4]=reshape(Float64[150.0; 1.64945 ],2,1) # 148.637
gt[:x5]=reshape(Float64[200.0; 1.77992 ],2,1) # 198.62
gt[:x6]=reshape(Float64[240.0; 2.20466 ],2,1) # 238.492
gt[:x7]=reshape(Float64[300.0; 2.14353 ],2,1) # 298.467
gt[:l1]=reshape(Float64[165.0; 1.17284 ],2,1) # 164.102


tree = prepBatchTree!(fg, drawpdf=false);

# list vertices in fg
@show xx,ll = ls(fg)

# do belief propagation inference over tree once
# using recursive single core approach (better stack trace for development)
# inferOverTreeR!(fg, tree)
inferOverTreeR!(fg, tree, N=N, dbg=true)
#
# test multi-processor solve (operational fast solving)
inferOverTree!(fg, tree)
# inferOverTree!(fg, tree, dbg=true)





#
