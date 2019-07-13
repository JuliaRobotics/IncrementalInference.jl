# addprocs(3)
using Test
# using Compat
# using IncrementalInference

@testset "out of module evalPotential..." begin
    include("TestModuleFunctions.jl")
end

include("testStateMachine.jl")

include("testCompareVariablesFactors.jl")

@testset "Ensure memory return is working properly..." begin
    include("typeReturnMemRef.jl")
end

include("basicGraphsOperations.jl")

include("testPartialFactors.jl")

@testset "basic Bayes tree construction" begin
    include("testBayesTreeiSAM2Example.jl")
end
@testset "Ensure converter types can be run from extending namespaces..." begin
    include("saveconvertertypes.jl")
end
@testset "packing converters work..." begin
    include("testgraphpackingconverters.jl")
end
include("testNLsolve.jl")

@testset "generic root finding by numeric solve of residual functions..." begin
    include("testNumericRootGenericRandomized.jl")
end
@testset "GenericWrapParam functors..." begin
    include("testCommonConvWrapper.jl")
end
@testset "with simple local constraint examples Odo, Obsv2..." begin
    include("testlocalconstraintexamples.jl")
end

include("testFactorMetadata.jl")

include("testExplicitMultihypo.jl")

include("testMultiHypo2Door.jl")

include("testMultithreaded.jl")

@testset "partial constraints..." begin
    include("testpartialconstraint.jl")
end
@testset "null hypothesis..." begin
    include("testnullhypothesis.jl")
end
@testset "standardized multihypothesis..." begin
    include("testmultihypothesisapi.jl")
end
@testset "with local Graphs.jl dictionary and arrays only (multicore)..." begin
    include("fourdoortest.jl")
end

@testset "saving to and loading from .jld2 file" begin

# DFG.savedfg(fg)
# DFG.loaddfg() # some convert problem on DFGVariable

savejld(fg, file="tempfg.jld2" )
@warn "not able to load a new DFG object yet"
# fgu = loadjld( file="tempfg.jld2" )
Base.rm("tempfg.jld2")

end

@warn "must return testExpandedJLD.jl to testing -- currently skipped since jld2 files cannot be loaded."
# include("testExpandedJLD.jl")












#
