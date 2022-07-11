
using IncrementalInference
using Statistics
using Test
using TensorCast

##

N=100
fg = initfg()

doors = [[0.0;],]
pd = manikde!(ContinuousScalar, doors; bw=[3.0;])
pd = resample(pd, N);
bws = getBW(pd)[:,1]
doors2 = getPoints(pd);

##

@testset "test evaluation of pose pose constraint" begin

##

v1 = addVariable!(fg,:x1, ContinuousScalar,N=N)
f1  = addFactor!(fg,[v1], Prior( pd )) #, samplefnc=getSample

# tem = 2.0*randn(1,N)+getVal(v1)+50.0
v2 = addVariable!(fg,:x2, ContinuousScalar, N=N)
odoc = LinearRelative(Normal(50.0,2.0)) # Odo(50.0*ones(1,1),2.0*ones(1,1),[1.0])
f2 = addFactor!(fg, [:x1; :x2], odoc ) #, samplefnc=getSample

# @test isInitialized(fg, :x1)

pts_ = approxConv(fg, :x1x2f1, :x2)
@cast pts[i,j] := pts_[j][i]
# pts = evalFactor(fg, f2, v2.label)
@show Statistics.mean(pts,dims=2)
@test norm(Statistics.mean(pts,dims=2)-[50.0]) < 15.0

tree = solveTree!(fg)

pts_ = getVal(fg, :x2)
@cast pts[i,j] := pts_[j][i]
@test norm(Statistics.mean(pts,dims=2)-[50.0]) < 15.0

##

end

# using RoMEPlotting
#
# plotKDE( getBelief(fg,:x2) )
# plotKDE(kde!(pts))









#
