using Base.Test

println("[TEST] out of module evalPotential...")
include("TestModuleFunctions.jl")
println("Success")

print("[TEST] Ensure memory return is working properly...")
include("typeReturnMemRef.jl")
println("Success")

print("[TEST] Ensure converter types can be run from extending namespaces...")
include("saveconvertertypes.jl")
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

println("[TEST] with simple local constraint examples Odo, Obsv2...")
include("testlocalconstraintexamples.jl")
println("Success")

println("[TEST] partial constraints...")
include("testpartialconstraint.jl")
println("Success")

println("[TEST] partial constraints...")
include("testnullhypothesis.jl")
println("Success")

println("[TEST] with local Graphs.jl dictionary and arrays only (multicore)...")
include("fourdoortest.jl")
println("Success")

println("[TEST] saving to and loading from .jld file")
savejld(fg, file="tempfg.jld" )
fgu = loadjld( file="tempfg.jld" )
Base.rm("tempfg.jld")
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


warn("plotMCMC needs ImageMagick on osx, not running test yet.")
# plotMCMC(tree, :x1, show=false)
# println("Success")














#
