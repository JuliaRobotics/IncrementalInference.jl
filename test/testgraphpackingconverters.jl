# packing and unpacking of Graph related types


using IncrementalInference, KernelDensityEstimate
using Base.Test

fg = emptyFactorGraph()

N=100

doors = [-100.0;0.0;100.0;300.0]'
pd = kde!(doors,[3.0])
pd = resample(pd,N);
bws = getBW(pd)[:,1]
doors2 = getPoints(pd);
v1 = addNode!(fg,:x1,doors,N=N)
f1  = addFactor!(fg,[v1],Obsv2( doors2, bws', [1.0]))

tem = 2.0*randn(1,N)+getVal(v1)+50.0
v2 = addNode!(fg,:x2, tem, N=N)
addFactor!(fg, [:x1; :x2], Odo([50.0]',[2.0]',[1.0]))


# using fourdoortest data
topack = getData(fg.g.vertices[4]) #fg.f[4].attributes["data"]
dd = convert(PackedFunctionNodeData{PackedOdo},topack)
upd = convert(FunctionNodeData{Odo}, dd)

@test topack.fnc.usrfnc!.Zij[1] == upd.fnc.usrfnc!.Zij[1]

# data structure conversion tests for protobuffing
println("Testing conversion to packed data structure and back")
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
