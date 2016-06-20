using IncrementalInference, KernelDensityEstimate
using Base.Test

include("fourdoortest.jl")
@test true


# using fourdoortest data
topack = fg.f[4].attributes["data"]
dd = convert(FunctionNodeData{PackedOdo},topack)
upd = convert(FunctionNodeData{Odo}, dd)
@test topack.fnc.Zij[1] == upd.fnc.Zij[1]

# data structure conversion tests for protobuffing
println("Testing conversion to packed data structure and back")
pd = convert(PackedVariableNodeData,fg.v[1].attributes["data"])
unpckd = convert(VariableNodeData, pd)

@test compare(fg.v[1].attributes["data"], unpckd)
@test fg.v[1].attributes["data"] == unpckd
println("Conversions and comparisons agree")



println("Run plot functions")
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
