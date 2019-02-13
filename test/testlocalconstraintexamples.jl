# test Odo, Obsv2

using IncrementalInference
using Statistics
using Test




global N=100
global fg = FactorGraph()

global doors = zeros(1,1)
global pd = kde!(doors,[3.0])
global pd = resample(pd,N);
global bws = getBW(pd)[:,1]
global doors2 = getPoints(pd);


@testset "test evaluation of pose pose constraint" begin

  global v1 = addVariable!(fg,:x1, ContinuousScalar,N=N)
  global f1  = addFactor!(fg,[v1], Prior( pd )) #, samplefnc=getSample

  # tem = 2.0*randn(1,N)+getVal(v1)+50.0
  global v2 = addVariable!(fg,:x2, ContinuousScalar, N=N)
  global odoc = LinearConditional(Normal(50.0,2.0)) # Odo(50.0*ones(1,1),2.0*ones(1,1),[1.0])
  global f2 = addFactor!(fg, [:x1; :x2], odoc ) #, samplefnc=getSample

  @test isInitialized(fg, :x1)

  global pts = evalFactor2(fg, f2, v2.index)
  @show Statistics.mean(pts,dims=2)
  @test norm(Statistics.mean(pts,dims=2)-[50.0]) < 15.0

  ensureAllInitialized!(fg)
  global tree = wipeBuildNewTree!(fg, drawpdf=false)
  inferOverTree!(fg, tree)

  @test norm(Statistics.mean(getVal(fg, :x2),dims=2)-[50.0]) < 15.0
end

# using RoMEPlotting
#
# plotKDE( getVertKDE(fg,:x2) )
# plotKDE(kde!(pts))









#
