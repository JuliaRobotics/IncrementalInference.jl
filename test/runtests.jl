# addprocs(3)

using Test
# using Compat
# using IncrementalInference

println("[TEST] out of module evalPotential...")
include("TestModuleFunctions.jl")
println("Success")

print("[TEST] Ensure memory return is working properly...")
include("typeReturnMemRef.jl")
println("Success")

println("[TEST] basic Bayes tree construction")
include("testBayesTreeiSAM2Example.jl")
println("Success")

print("[TEST] Ensure converter types can be run from extending namespaces...")
include("saveconvertertypes.jl")
println("Success")

println("[TEST] packing converters work...")
include("testgraphpackingconverters.jl")
println("Success")

include("testNLsolve.jl")

println("[TEST] generic root finding by numeric solve of residual functions...")
include("testNumericRootGenericRandomized.jl")
println("Success")

println("[TEST] GenericWrapParam functors...")
include("testCommonConvWrapper.jl")
println("Success")

println("[TEST] with simple local constraint examples Odo, Obsv2...")
include("testlocalconstraintexamples.jl")
println("Success")

include("testFactorMetadata.jl")

include("testExplicitMultihypo.jl")

include("testMultiHypo2Door.jl")

include("testMultithreaded.jl")

println("[TEST] partial constraints...")
include("testpartialconstraint.jl")
println("Success")

println("[TEST] null hypothesis...")
include("testnullhypothesis.jl")
println("Success")

println("[TEST] standardized multihypothesis...")
include("testmultihypothesisapi.jl")
println("Success")

println("[TEST] with local Graphs.jl dictionary and arrays only (multicore)...")
include("fourdoortest.jl")
println("Success")

println("[TEST] saving to and loading from .jld2 file")
savejld(fg, file="tempfg.jld2" )
fgu = loadjld( file="tempfg.jld2" )
Base.rm("tempfg.jld2")
println("Success")

include("testExpandedJLD.jl")














#
