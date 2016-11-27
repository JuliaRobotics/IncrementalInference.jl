using IncrementalInference, TransformUtils
using Base.Test

using KernelDensityEstimate

N = 200
fg = emptyFactorGraph()

initCov = eye(6)
[initCov[i,i] = 0.01 for i in 4:6];
odoCov = deepcopy(initCov)


println("Adding PriorPose3 to graph...")
v1 = addNode!(fg,"x1",  0.1*randn(6,N),  N=N)
initPosePrior = PriorPose3(SE3(0), initCov)
f1  = addFactor!(fg,[v1], initPosePrior)

println("Ensure vertex initialized properly")
# start with to tight an initialization
muX1 = Base.mean(getVal(fg,"x1"),2)
stdX1 = Base.std(getVal(fg,"x1"),2)
@test sum(map(Int,abs(muX1) .< 0.1)) == 6
@test sum(map(Int, 0.05 .< stdX1 .< 0.15)) == 6


println("Testing PriorPose3 evaluation...")
priorpts = evalFactor2(fg, fg.g.vertices[2], 1)
means = Base.mean(priorpts,2)
@test sum(map(Int,abs(means[1:3]) .> 0.5)) == 0
@test sum(map(Int,abs(means[4:6]) .> 0.05)) == 0


println("Adding Pose3Pose3 to graph...")
odo = SE3([10;0;0], Quaternion(0))
pts0X2 = projectParticles(getVal(fg,"x1"), odo, odoCov)
odoconstr = Pose3Pose3(odo, odoCov)
v2 = addNode!(fg,"x2",  pts0X2, N=N)
addFactor!(fg,[v1;v2],odoconstr)


println("Testing Pose3Pose3 evaluation...")
X1pts = evalFactor2(fg, fg.g.vertices[4], 1)
X2pts = evalFactor2(fg, fg.g.vertices[4], 3)
X2ptsMean = Base.mean(X2pts,2)
X1ptsMean = Base.mean(X1pts,2)
@test  sum(map(Int, abs(X1ptsMean) .< 0.5 )) == 6
@test  sum(map(Int, abs(X2ptsMean - [10.0;0;0;0;0;0]) .< 0.5 )) == 6


println("Construct Bayes tree and perform inference...")
tree = prepBatchTree!(fg);
inferOverTree!(fg, tree)

println("Ensure basic parameters on x1,x2 after inference...")
# check mean and covariances after one up and down pass over the tree
muX1 = Base.mean(getVal(fg,"x1"),2)
stdX1 = Base.std(getVal(fg,"x1"),2)
@test sum(map(Int,abs(muX1[1:3]) .< 1.0)) == 3
@test sum(map(Int,abs(muX1[4:6]) .< 0.1)) == 3
@test sum(map(Int, 0.5 .< stdX1[1:3] .< 1.5)) == 3
@test sum(map(Int, 0.025 .< stdX1[4:6] .< 0.25)) == 3
muX2 = Base.mean(getVal(fg,"x2"),2)
stdX2 = Base.std(getVal(fg,"x2"),2)
@test sum(map(Int, abs(muX2[1:3]-[10.0;0;0]) .< 1.0)) == 3
@test sum(map(Int, abs(muX2[4:6]) .< 0.1)) == 3
@test sum(map(Int, 1.0 .< stdX2[1:3] .< 2.0)) == 3
@test sum(map(Int, 0.05 .< stdX2[4:6] .< 0.25)) == 3

# println("Plot marginals to see what is happening")
# plotKDE(marginal(getVertKDE(fg,"x1"),[1]))
# plotKDE(marginal(getVertKDE(fg,"x2"),[1]))

println("Modify factor graph, adding Pose3Pose3's X3,X4 to graph...")
#X3, yaw 90 degrees
odo2 = SE3([0;0;0], AngleAxis(pi/2,[0;0;1]))
odoconstr2 = Pose3Pose3(odo2, odoCov)
pts0x3 = projectParticles(getVal(fg,"x2"), odo2, odoCov) # init new pose location
v3 = addNode!(fg,"x3", pts0x3, N=N)
addFactor!(fg,[v2;v3],odoconstr2)

# X4, drive forward 10m
odo3 = SE3([10;0;0], SO3(0))
odoconstr3 = Pose3Pose3(odo3, odoCov)
pts0X4 = projectParticles(getVal(fg,"x3"),odo3, odoCov)
v4 = addNode!(fg,"x4", pts0X4, N=N)
addFactor!(fg,[v3;v4],odoconstr3)

# println("Plot marginals to see what is happening")
# plotKDE(marginal(getVertKDE(fg,"x3"),[1]))
# plotKDE(marginal(getVertKDE(fg,"x4"),[1]))

println("Reconstruct Bayes tree and perform inference...")
tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)




#
