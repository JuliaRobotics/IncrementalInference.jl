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
    @testset "Temporary Debug Group (incl. frequent numerical issues)" begin
    include("testBasicGraphs.jl") # NUMERICAL
    include("testEuclidDistance.jl") # test_broken
    include("testSpecialSampler.jl") # FIX, bounds error
    include("testDERelative.jl") # FIX BoundsError Ln332 cf._legacyParams[k][i], [100] of 1..99
    include("testHeatmapGridDensity.jl") # numerical test broken
    include("testMultiHypo3Door.jl") # FIX numerical
    include("testExpXstroke.jl") # FIX, init bw issue, addLikelihoodsDifferentialCHILD! LN327
    include("testCircular.jl") # FIX
    include("testFluxModelsDistribution.jl")
    include("testMixturePrior.jl") # FIX, serde structutils issue with BinarTruckFixedDepth?
    include("testMixtureLinearConditional.jl") # FIX, use HomotopyDensity as replacement for Mixture
    include("testMultihypoAndChain.jl") # numerical issues

    @error "See new DFG v0.29 JSON serde for variables and factors, old IIF serde tests currently disabled"
    if false
      # include("TestModuleFunctions.jl")
      include("testCompareVariablesFactors.jl")
      include("saveconvertertypes.jl")
      include("testgraphpackingconverters.jl")
      include("testSaveLoadDFG.jl")
      include("testPackingMixtures.jl")
    end
  end
end

if TEST_GROUP in ["all", "basic_functional_group"]
@testset "Basic Functional Group" begin

# start as basic as possible and build from there
include("typeReturnMemRef.jl")
include("testDistributionsGeneric.jl")
include("basicGraphsOperations.jl")

#FIXME fails on MetaBayesTree
include("testTreeSaveLoad.jl")

# test convolution functions
include("testApproxConv.jl")
include("testBasicForwardConvolve.jl")

include("testDefaultDeconv.jl")
include("priorusetest.jl") # many skips

include("testCliqSolveDbgUtils.jl")
include("testUseMsgLikelihoods.jl")

@test_broken error("testSphereMani.jl broken") # include("testSphereMani.jl") # FIXME
include("testBasicManifolds.jl")
include("testSpecialOrthogonalMani.jl")

# gradient / jacobian tests
#include("manifolds/manifolddiff.jl")
#include("manifolds/factordiff.jl")
@error "Gradient tests must be updated and restored for new ccw.varValsAll[]"
#include("testGradientUtils.jl")
#include("testFactorGradients.jl")

include("testJunctionTreeConstruction.jl")
include("testBayesTreeiSAM2Example.jl")
include("testTreeFunctions.jl")

include("testCommonConvWrapper.jl") # skipped test with just ::Float point type 

include("testStateMachine.jl")
include("testBasicCSM.jl")
include("testCliqueFactors.jl")
include("testCcolamdOrdering.jl")
include("testCliqueTreesOrderings.jl")
include("testJointEnforcement.jl")
include("testHasPriors913.jl")
include("testInitVariableOrder.jl")
include("testTreeMessageUtils.jl")
include("testCSMMonitor.jl")
include("testBasicRecycling.jl")
include("testSkipUpDown.jl")
include("testlocalconstraintexamples.jl")

include("testManualInit.jl")
include("testBasicTreeInit.jl")
include("testSolveOrphanedFG.jl")
include("testSolveKey.jl")

include("testPartialFactors.jl")
include("testPartialPrior.jl")

## WORK IN PROGRESS
include("testSpecialEuclidean2Mani.jl") # TBD
include("testpartialconstraint.jl") # FIX
include("testPartialNH.jl") # FIX
include("testBasicParametric.jl") # FIX

end
end

if TEST_GROUP in ["all", "test_cases_group"]
@testset "Test Cases Group" begin

include("testnullhypothesis.jl") 
include("testVariousNSolveSize.jl")
include("testExplicitMultihypo.jl")
include("TestCSMMultihypo.jl")
include("testCalcFactorHypos.jl")
include("testMultimodal1D.jl") # skipped a test
include("testMultithreaded.jl")
include("testmultihypothesisapi.jl")
include("fourdoortest.jl")
include("testAnalysisTools.jl")


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
