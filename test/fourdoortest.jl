using IncrementalInference

@testset "test fourdoor early example" begin

global N=100
global fg = initfg()
global doors = reshape(Float64[-100.0;0.0;100.0;300.0],1,4)
global pd = kde!(doors,[3.0])
global pd = resample(pd,N);
global bws = getBW(pd)[:,1]
global doors2 = getPoints(pd);


addVariable!(fg,:x0,ContinuousScalar,N=N)
addFactor!(fg,[:x0], Prior( pd ) )

# tem = 2.0*randn(1,N)+getVal(v1)+50.0
addVariable!(fg,:x2, ContinuousScalar, N=N)
addFactor!(fg, [:x0; :x2], LinearConditional(Normal(50.0,2.0)) )
# addFactor!(fg, [v1;v2], Odo(50.0*ones(1,1),[2.0]',[1.0]))


# monocular sighting would look something like
#addFactor!(fg, Mono, [:x3,:l1], [14.0], [1.0], [1.0])
#addFactor!(fg, Mono, [:x4,:l1], [11.0], [1.0], [1.0])

addVariable!(fg,:x3,ContinuousScalar, N=N)
addFactor!(fg,[:x2;:x3], LinearConditional( Normal(50.0,4.0)) )
addFactor!(fg,[:x3], Prior( pd ))

addVariable!(fg,:x4,ContinuousScalar, N=N)
addFactor!(fg,[:x3;:x4], LinearConditional( Normal(50.0,2.0)) )


addVariable!(fg, :l1, ContinuousScalar, N=N)
addFactor!(fg, [:x3,:l1], Ranged([64.0],[0.5],[1.0]))
addFactor!(fg, [:x4,:l1], Ranged([16.0],[0.5],[1.0]))



addVariable!(fg,:x5,ContinuousScalar, N=N)
addFactor!(fg,[:x4;:x5], LinearConditional( Normal(50.0,2.0)) )


addVariable!(fg,:x6,ContinuousScalar, N=N)
addFactor!(fg,[:x5;:x6], LinearConditional( Normal(40.0,1.20)) )


addVariable!(fg,:x7,ContinuousScalar, N=N)
addFactor!(fg,[:x6;:x7], LinearConditional( Normal(60.0,2.0)) )

# ensureAllInitialized!(fg)

global mlc = MixturePrior(Normal.(doors[1,:], bws[1]), 0.25*ones(4))

# getSample(mlc)

addFactor!(fg,[:x7], mlc )



# HMM computed ground truth, extended for 7 poses with landmark
global gt = Dict{Symbol, Array{Float64,2}}()
gt[:x0]=reshape(Float64[0.0;1.97304 ],2,1) # -0.0342366
gt[:x2]=reshape(Float64[50.0; 2.83153 ],2,1) # 49.8797
gt[:x3]=reshape(Float64[100.0; 1.65557 ],2,1) # 99.8351
gt[:x4]=reshape(Float64[150.0; 1.64945 ],2,1) # 148.637
gt[:x5]=reshape(Float64[200.0; 1.77992 ],2,1) # 198.62
gt[:x6]=reshape(Float64[240.0; 2.20466 ],2,1) # 238.492
gt[:x7]=reshape(Float64[300.0; 2.14353 ],2,1) # 298.467
gt[:l1]=reshape(Float64[165.0; 1.17284 ],2,1) # 164.102


tree, smt, hist = solveTree!(fg)

# global tree = prepBatchTree!(fg, drawpdf=false);
#
# # list vertices in fg
# @show xx,ll = IncrementalInference.ls(fg)
#
# # do belief propagation inference over tree once
# # using recursive single core approach (better stack trace for development)
# # inferOverTreeR!(fg, tree)
# inferOverTree!(fg, tree, N=N, dbg=true, treeinit=true)
#
#  #
# # test multi-processor solve (operational fast solving)
#  inferOverTree!(fg, tree)
# # inferOverTree!(fg, tree, dbg=true)
#


end

#
