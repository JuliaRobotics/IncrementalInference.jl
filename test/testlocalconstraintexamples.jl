# test Odo, Obsv2

using IncrementalInference, KernelDensityEstimate
using Base.Test



fg = emptyFactorGraph()

N=100

doors = zeros(1,1)
pd = kde!(doors,[3.0])
pd = resample(pd,N);
bws = getBW(pd)[:,1]
doors2 = getPoints(pd);

@testset "test evaluation of pose pose constraint" begin

  v1 = addNode!(fg,:x1, ContinuousScalar,N=N)
  f1  = addFactor!(fg,[v1],Obsv2( doors2, bws', [1.0])) #, samplefnc=getSample

  # tem = 2.0*randn(1,N)+getVal(v1)+50.0
  v2 = addNode!(fg,:x2, ContinuousScalar, N=N)
  odoc = Odo(50.0*ones(1,1),2.0*ones(1,1),[1.0])
  f2 = addFactor!(fg, [:x1; :x2], odoc ) #, samplefnc=getSample

  @test isInitialized(fg, :x1)

  pts = evalFactor2(fg, f2, v2.index)
  @show Base.mean(pts,2)
  @test norm(Base.mean(pts,2)-[50.0]) < 15.0

  ensureAllInitialized!(fg)
  tree = wipeBuildNewTree!(fg, drawpdf=false)
  inferOverTree!(fg, tree)

  @test norm(Base.mean(getVal(fg, :x2),2)-[50.0]) < 15.0
end

# using RoMEPlotting
#
# plotKDE( getVertKDE(fg,:x2) )
# plotKDE(kde!(pts))













#
