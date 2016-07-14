using Base.Test

addprocs(3)

println("[TEST] with local Graphs.jl dictionary and arrays only (multicore)...")
include("fourdoortest.jl")
println("Success")
@test true


println("[TEST] plot functions...")
using Gadfly
# draw all beliefs
DOYTICKS = false
xx,ll = ls(fg)
msgPlots = drawHorBeliefsList(fg, xx, gt=gt,nhor=2);
evalstr = ""
for i in 1:length(msgPlots)
    evalstr = string(evalstr, ",msgPlots[$(i)]")
end
pl = eval(parse(string("vstack(",evalstr[2:end],")")));
@test true

print("[TEST] Ensure memory return is working properly...")
include("typeReturnMemRef.jl")
println("Success")


println("[TEST] packing converters work...")
# using fourdoortest data
topack = fg.f[4].attributes["data"]
dd = convert(FunctionNodeData{PackedOdo},topack)
upd = convert(FunctionNodeData{Odo}, dd)
@test topack.fnc.Zij[1] == upd.fnc.Zij[1]

# data structure conversion tests for protobuffing
println("Testing conversion to packed data structure and back")
dat = IncrementalInference.dlapi.getvertex(fg,1).attributes["data"]
pd = convert(PackedVariableNodeData, dat) #fg.v[1].attributes["data"]
# pd = convert(PackedVariableNodeData,fg.v[1].attributes["data"])
unpckd = convert(VariableNodeData, pd)

@test compare(IncrementalInference.dlapi.getvertex(fg,1).attributes["data"], unpckd)
@test IncrementalInference.dlapi.getvertex(fg,1).attributes["data"] == unpckd
println("Conversions and comparisons agree")



println("[TEST] with CloudGraphs data layer (multicore)...")
include("fourdoortestcloudgraph.jl")
println("Success")
@test true
