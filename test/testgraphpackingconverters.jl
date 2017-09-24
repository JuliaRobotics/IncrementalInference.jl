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
v1 = addNode!(fg,:x1,doors,N=N)
f1  = addFactor!(fg,[v1],Obsv2( doors2, bws2, [1.0]))

tem = 2.0*randn(1,N)+getVal(v1)+50.0
v2 = addNode!(fg,:x2, tem, N=N)
f2 = addFactor!(fg, [:x1, :x2], Odo(50.0*ones(1,1),2.0*ones(1,1),[1.0]))


println("Testing conversion to packed function node data structure and back")

topack = getData(f1) #fg.f[4].attributes["data"]
dd = convert(PackedFunctionNodeData{PackedObsv2},topack)
upd = convert(FunctionNodeData{GenericWrapParam{Obsv2}}, dd)

@test compare(topack, upd)

topack = getData(f2) #fg.f[4].attributes["data"]
dd = convert(PackedFunctionNodeData{PackedOdo},topack)
upd = convert(FunctionNodeData{GenericWrapParam{Odo}}, dd)

@test compare(topack, upd)


packedv4data = FNDencode(IncrementalInference.PackedFunctionNodeData{PackedOdo}, getData(f2))
upv4data = FNDdecode(IncrementalInference.FunctionNodeData{GenericWrapParam{Odo}}, packedv4data)

@test compare(getData(f2), upv4data)

# data structure conversion tests for protobuffing
println("Testing conversion to packed variable node data structure and back")
dat = getData(fg,1)
pd = convert(PackedVariableNodeData, dat) #fg.v[1].attributes["data"]
# pd = convert(PackedVariableNodeData,fg.v[1].attributes["data"])
unpckd = convert(VariableNodeData, pd)

@test IncrementalInference.compare(getData(fg,1), unpckd)
@test getData(fg,1) == unpckd


packedv1data = VNDencoder(IncrementalInference.PackedVariableNodeData, getData(v1))
upv1data = VNDdecoder(IncrementalInference.VariableNodeData, packedv1data)

@test compare(getData(v1), upv1data)

println("Conversions and comparisons agree")






















#
