using IncrementalInference, KernelDensityEstimate
using Base.Test

include("fourdoortest.jl")
@test true


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
