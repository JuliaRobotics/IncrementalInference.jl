# addprocs(3)
using Test
# using Compat
# using IncrementalInference

@testset "Parametric Tests" begin
    include("testBasicParametric.jl")
end

@testset "out of module evalPotential..." begin
    include("TestModuleFunctions.jl")
end

include("testStateMachine.jl")

include("testCompareVariablesFactors.jl")

@testset "Ensure memory return is working properly..." begin
    include("typeReturnMemRef.jl")
end

include("basicGraphsOperations.jl")

include("testMixturePrior.jl")

include("testPartialFactors.jl")

@testset "basic Bayes tree construction" begin
    include("testBayesTreeiSAM2Example.jl")
end

#FIXME fails on MetaBayesTree
include("testTreeSaveLoad.jl")

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

include("testBasicForwardConvolve.jl")

include("testFactorMetadata.jl")

include("testJunctionTreeConstruction.jl")

include("testBasicCSM.jl")

include("testCliqueFactors.jl")

include("testCcolamdOrdering.jl")

include("testBasicGraphs.jl")

@testset "with simple local constraint examples Odo, Obsv2..." begin
    # old names should be removed, like Odo, Obsv2
    include("testlocalconstraintexamples.jl")
end

include("testSolveOrphanedFG.jl")

include("testSolveSetPPE.jl")

# include("priorusetest.jl")

include("testVariousNSolveSize.jl")

include("testExplicitMultihypo.jl")

include("TestCSMMultihypo.jl")

include("testMultiHypo2Door.jl")

include("testMultimodal1D.jl")

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

include("testAnalysisTools.jl")

@testset "saving to and loading from FileDFG" begin
    saveFolder = "/tmp/dfg_test"
    # VERSION above 1.0.x hack required since Julia 1.0 does not seem to have function `splitpath`
    saveDFG(fg, saveFolder, compress= VERSION < v"1.1" ? :none : :gzip)
    retDFG = LightDFG{SolverParams}(params=SolverParams())
    if v"1.1" <= VERSION
      retDFG = loadDFG(saveFolder, IncrementalInference, retDFG)
      @test symdiff(ls(fg), ls(retDFG)) == []
      @test symdiff(lsf(fg), lsf(retDFG)) == []
    end
end


# dont run test on ARM, as per issue #527
if Base.Sys.ARCH in [:x86_64;]
  include("testTexTreeIllustration.jl")
end










#
