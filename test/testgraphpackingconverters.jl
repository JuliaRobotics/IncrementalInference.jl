# packing and unpacking of Graph related types


using KernelDensityEstimate
using IncrementalInference

using Base.Test

fg = emptyFactorGraph()

N=100

doors = zeros(1,4)
doors[1,1:4] = [-100.0;0.0;100.0;300.0]
pd = kde!(doors,[3.0])
pd = resample(pd,N);
bws = getBW(pd)[:,1]
bws2 = Array{Float64}(length(bws),1)
bws2[:,1] = bws[:]
doors2 = getPoints(pd);
v1 = addNode!(fg,:x1, ContinuousScalar,N=N)
f1  = addFactor!(fg,[v1],Obsv2( doors2, bws2, [1.0]))

tem = 2.0*randn(1,N)+getVal(v1)+50.0
v2 = addNode!(fg,:x2, ContinuousScalar, N=N)
f2 = addFactor!(fg, [:x1, :x2], Odo(50.0*ones(1,1),2.0*ones(1,1),[1.0]))


@testset "Testing conversion to packed function node data structure and back" begin
  topack = getData(f1)
  dd = convert(PackedFunctionNodeData{PackedObsv2},topack)
  upd = convert(FunctionNodeData{CommonConvWrapper{Obsv2}}, dd)

  @test compare(topack, upd)

  topack = getData(f2) #fg.f[4].attributes["data"]
  dd =  convert(IncrementalInference.PackedFunctionNodeData{PackedOdo},topack)
  upd = convert(IncrementalInference.FunctionNodeData{CommonConvWrapper{Odo}}, dd)

  @test compare(topack, upd)

  # packedv4data = FNDencode(IncrementalInference.PackedFunctionNodeData{PackedOdo}, getData(f2))
  # upv4data = FNDdecode(IncrementalInference.FunctionNodeData{CommonConvWrapper{Odo}}, packedv4data)
  # @test compare(getData(f2), upv4data)
end

@testset "Testing conversion to packed variable node data structure and back" begin
  dat = getData(fg,:x1) # 1
  pd = convert(IncrementalInference.PackedVariableNodeData, dat)
  unpckd = convert(IncrementalInference.VariableNodeData, pd)

  @test IncrementalInference.compare(getData(fg, :x1), unpckd)
  @test getData(fg,:x1) == unpckd

  # packedv1data = VNDencoder(IncrementalInference.PackedVariableNodeData, getData(v1))
  # upv1data = VNDdecoder(IncrementalInference.VariableNodeData, packedv1data)
  # @test compare(getData(v1), upv1data)
end






















#
