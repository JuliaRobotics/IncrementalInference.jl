##

using Test

# TODO remove, forcing conflict to use LieGroups
using LieGroups: TranslationGroup
using DistributedFactorGraphs
DFG.@usingDFG true

TEST_GROUP = get(ENV, "IIF_TEST_GROUP", "all")

##

@testset "IncrementalInference Tests" begin
  # temporarily moved to start (for debugging)
  if TEST_GROUP in ["all", "tmp_debug_group"]
    @testset "Temporary Debug Group" begin
    include("testSpecialOrthogonalMani.jl")
    include("testMultiHypo3Door.jl") # FIX numerical
    include("priorusetest.jl") # MANY SKIPS
  end
end

if TEST_GROUP in ["all", "basic_functional_group"]
@testset "Basic Functional Group" begin

# start as basic as possible and build from there
include("typeReturnMemRef.jl")
include("testDistributionsGeneric.jl")
include("basicGraphsOperations.jl")

@error "See new DFG v0.29 JSON serde for variables and factors, old IIF serde tests currently disabled"
if false
  # include("TestModuleFunctions.jl")
  include("testCompareVariablesFactors.jl")
  include("saveconvertertypes.jl")
  include("testgraphpackingconverters.jl")
  include("testSaveLoadDFG.jl")
  include("testPackingMixtures.jl")
end
#FIXME fails on MetaBayesTree
include("testTreeSaveLoad.jl")

# test convolution functions
include("testApproxConv.jl")
include("testBasicForwardConvolve.jl")

include("testDefaultDeconv.jl")
include("testCliqSolveDbgUtils.jl")
include("testUseMsgLikelihoods.jl")

include("testEuclidDistance.jl") # FIX
@test_broken error("testSphereMani.jl broken") # include("testSphereMani.jl") # FIXME
include("testBasicManifolds.jl")
include("testSpecialEuclidean2Mani.jl") # TBD
# gradient / jacobian tests
#include("manifolds/manifolddiff.jl")
#include("manifolds/factordiff.jl")
@error "Gradient tests must be updated and restored for new ccw.varValsAll[]"
#include("testGradientUtils.jl")
#include("testFactorGradients.jl")

include("testJunctionTreeConstruction.jl")
include("testBayesTreeiSAM2Example.jl")
include("testTreeFunctions.jl")

include("testCommonConvWrapper.jl") # FIX
include("testSpecialSampler.jl") # TODO, rename, refine
include("testHeatmapGridDensity.jl") # FIX

include("testStateMachine.jl")
include("testBasicCSM.jl")
include("testCliqueFactors.jl")
include("testCcolamdOrdering.jl")
include("testCliqueTreesOrderings.jl")
include("testBasicGraphs.jl") # NUMERICAL
include("testJointEnforcement.jl") # FIX
include("testHasPriors913.jl") # FIX something in solve
include("testInitVariableOrder.jl")
include("testTreeMessageUtils.jl")
include("testCSMMonitor.jl")
include("testExpXstroke.jl") # FIX
include("testBasicRecycling.jl") # FIX
include("testSkipUpDown.jl") # FIX
include("testlocalconstraintexamples.jl")
include("testManualInit.jl") # FIX
include("testBasicTreeInit.jl") # FIX
include("testSolveOrphanedFG.jl")
include("testSolveKey.jl")

include("testPartialFactors.jl")
include("testPartialPrior.jl") # FIX
include("testpartialconstraint.jl") # FIX
include("testPartialNH.jl") # FIX
include("testMixturePrior.jl") # FIX

include("testDERelative.jl") # FIX BoundsError Ln332 cf._legacyParams[k][i], [100] of 1..99

end
end

if TEST_GROUP in ["all", "test_cases_group"]
@testset "Test Cases Group" begin
include("testnullhypothesis.jl") 
include("testVariousNSolveSize.jl")
include("testExplicitMultihypo.jl")
include("TestCSMMultihypo.jl")
include("testCalcFactorHypos.jl")
include("testMultimodal1D.jl") # FIX
include("testMultihypoAndChain.jl") # FIX numerical
include("testMultithreaded.jl")
include("testmultihypothesisapi.jl") # FIX
include("fourdoortest.jl")
include("testCircular.jl") # FIX
include("testMixtureLinearConditional.jl") # FIX
include("testAnalysisTools.jl")
if false
  include("testFluxModelsDistribution.jl")
else
  # @error "Skipped testFluxModelsDistribution.jl"
  @test_skip("Skipped testFluxModelsDistribution.jl")
end


include("testBasicParametric.jl") # FIX
# include("testMixtureParametric.jl") #FIXME parametric mixtures #1787

# dont run test on ARM, as per issue #527
if Base.Sys.ARCH in [:x86_64;]
  include("testTexTreeIllustration.jl")
end

# include("testMultiprocess.jl")
include("testDeadReckoningTether.jl")
end
end
end

#
