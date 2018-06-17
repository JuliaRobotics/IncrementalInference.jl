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
v1 = addNode!(fg,:x1, ContinuousScalar,N=N)
f1  = addFactor!(fg,[v1],Obsv2( doors2, bws', [1.0])) #, samplefnc=getSample

tem = 2.0*randn(1,N)+getVal(v1)+50.0
v2 = addNode!(fg,:x2, ContinuousScalar, N=N)
odoc = Odo(50.0*ones(1,1),2.0*ones(1,1),[1.0])
f2 = addFactor!(fg, [:x1; :x2], odoc ) #, samplefnc=getSample



@testset "test evaluation of pose pose constraint" begin
  pts = evalFactor2(fg, f2, v2.index)
  @test norm(Base.mean(pts,2)-[50.0]) < 15.0


  tree = wipeBuildNewTree!(fg, drawpdf=false)

  inferOverTree!(fg, tree)
  
  @test norm(Base.mean(getVal(fg, :x2),2)-[50.0]) < 10.0
end

# plotKDE( getVertKDE(fg,:x2) )














#
