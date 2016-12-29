using IncrementalInference, KernelDensityEstimate

fg = emptyFactorGraph()

N=100

doors = [-100.0;0.0;100.0;300.0]'
pd = kde!(doors,[3.0])
pd = resample(pd,N);
bws = getBW(pd)[:,1]
doors2 = getPoints(pd);
v1 = addNode!(fg,:x1,doors,N=N)
f1  = addFactor!(fg,[v1],Obsv2( doors2, bws', [1.0]), samplefnc=getSample)

tem = 2.0*randn(1,N)+getVal(v1)+50.0
v2 = addNode!(fg,:x2, tem, N=N)
addFactor!(fg, [v1;v2], Odo([50.0]',[2.0]',[1.0]), samplefnc=getSample)

# monocular sighting would look something like
#addFactor!(fg, Mono, [:x3,:l1], [14.0], [1.0], [1.0])
#addFactor!(fg, Mono, [:x4,:l1], [11.0], [1.0], [1.0])

v3=addNode!(fg,:x3,4.0*randn(1,N)+getVal(v2)+50.0, N=N)
addFactor!(fg,[v2;v3],Odo([50.0]',[4.0]',[1.0]))
f2 = addFactor!(fg,[v3], Obsv2(doors2, bws', [1.0]))

v4=addNode!(fg,:x4,2.0*randn(1,N)+getVal(v3)+50.0, N=N)
addFactor!(fg,[v3;v4],Odo([50.0]',[2.0]',[1.0]))


    l1=addNode!(fg, :l1, 0.5*randn(1,N)+getVal(v3)+64.0, N=N)
    addFactor!(fg, [v3,l1], Ranged([64.0],[0.5],[1.0]))
    addFactor!(fg, [v4,l1], Ranged([16.0],[0.5],[1.0]))



v5=addNode!(fg,:x5,2.0*randn(1,N)+getVal(v4)+50.0, N=N)
addFactor!(fg,[v4;v5],Odo([50.0]',[2.0]',[1.0]))


v6=addNode!(fg,:x6,1.25*randn(1,N)+getVal(v5)+40.0, N=N)
addFactor!(fg,[v5;v6],Odo([40.0]',[1.25]',[1.0]))


v7=addNode!(fg,:x7,2.0*randn(1,N)+getVal(v6) +60.0, N=N)
addFactor!(fg,[v6;v7],Odo([60.0]',[2.0]',[1.0]))

f3 = addFactor!(fg,[v7], Obsv2(doors, bws', [1.0]))


# HMM computed ground truth, extended for 7 poses with landmark
gt = Dict{Symbol, Array{Float64,2}}()
gt[:x1]=([0.0;1.97304 ]')' # -0.0342366
gt[:x2]=([50.0; 2.83153 ]')' # 49.8797
gt[:x3]=([100.0; 1.65557 ]')' # 99.8351
gt[:x4]=([150.0; 1.64945 ]')' # 148.637
gt[:x5]=([200.0; 1.77992 ]')' # 198.62
gt[:x6]=([240.0; 2.20466 ]')' # 238.492
gt[:x7]=([300.0; 2.14353 ]')' # 298.467
gt[:l1]=([165.0; 1.17284 ]')' # 164.102


tree = prepBatchTree!(fg, drawpdf=false);

# list vertices in fg
@show xx,ll = ls(fg)

# do belief propagation inference over tree once
inferOverTreeR!(fg, tree)

inferOverTree!(fg, tree)
println("Inference finished")
