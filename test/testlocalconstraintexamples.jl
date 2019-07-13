# test Odo, Obsv2

using IncrementalInference
using Statistics
using Test




N=100
fg = initfg()

doors = zeros(1,1)
pd = kde!(doors,[3.0])
pd = resample(pd,N);
bws = getBW(pd)[:,1]
doors2 = getPoints(pd);


@testset "test evaluation of pose pose constraint" begin

v1 = addVariable!(fg,:x1, ContinuousScalar,N=N)
f1  = addFactor!(fg,[v1], Prior( pd )) #, samplefnc=getSample

# tem = 2.0*randn(1,N)+getVal(v1)+50.0
v2 = addVariable!(fg,:x2, ContinuousScalar, N=N)
odoc = LinearConditional(Normal(50.0,2.0)) # Odo(50.0*ones(1,1),2.0*ones(1,1),[1.0])
f2 = addFactor!(fg, [:x1; :x2], odoc ) #, samplefnc=getSample

# @test isInitialized(fg, :x1)

pts = approxConv(fg, :x1x2f1, :x2)
# pts = evalFactor2(fg, f2, v2.label)
@show Statistics.mean(pts,dims=2)
@test norm(Statistics.mean(pts,dims=2)-[50.0]) < 15.0

tree, smt, hist = solveTree!(fg)
# ensureAllInitialized!(fg)
# tree = wipeBuildNewTree!(fg, drawpdf=false)
# inferOverTree!(fg, tree)

@test norm(Statistics.mean(getVal(fg, :x2),dims=2)-[50.0]) < 15.0

end

# using RoMEPlotting
#
# plotKDE( getVertKDE(fg,:x2) )
# plotKDE(kde!(pts))









#
