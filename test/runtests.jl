using Base.Test

println("[TEST] out of module evalPotential...")
include("TestModuleFunctions.jl")
println("Success")

print("[TEST] Ensure memory return is working properly...")
include("typeReturnMemRef.jl")
println("Success")

addprocs(3)

println("[TEST] packing converters work...")
include("testgraphpackingconverters.jl")
println("Success")

println("[TEST] generic root finding by numeric solve of residual functions...")
include("testNumericRootGenericRandomized.jl")
println("Success")

println("[TEST] GenericWrapParam functors...")
include("testGenericWrapParam.jl")
println("Success")

println("[TEST] with local Graphs.jl dictionary and arrays only (multicore)...")
include("fourdoortest.jl")
println("Success")


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
println("Success")
















#
