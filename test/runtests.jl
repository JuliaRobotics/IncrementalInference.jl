using Test

TEST_GROUP = get(ENV, "IIF_TEST_GROUP", "all")

# temporarily moved to start (for debugging)
#...
if TEST_GROUP in ["all", "tmp_debug_group"]
include("testMultiHypo3Door.jl")
include("priorusetest.jl")
end

if TEST_GROUP in ["all", "basic_functional_group"]
# more frequent stochasic failures from numerics
include("testSpecialEuclidean2Mani.jl")
include("testEuclidDistance.jl")

# regular testing
include("testSphereMani.jl")
include("testSpecialOrthogonalMani.jl")
include("testBasicManifolds.jl")

# start as basic as possible and build from there
include("typeReturnMemRef.jl")
include("testDistributionsGeneric.jl")
include("testHeatmapGridDensity.jl")
include("testCliqSolveDbgUtils.jl")
include("basicGraphsOperations.jl")

include("TestModuleFunctions.jl")
include("testCompareVariablesFactors.jl")
include("saveconvertertypes.jl")
include("testgraphpackingconverters.jl")
include("testSaveLoadDFG.jl")

include("testPackingMixtures.jl")

include("testJunctionTreeConstruction.jl")
include("testBayesTreeiSAM2Example.jl")
include("testTreeFunctions.jl")

#FIXME fails on MetaBayesTree
include("testTreeSaveLoad.jl")

include("testGradientUtils.jl")
include("testFactorGradients.jl")
include("testSpecialSampler.jl") # TODO, rename, refine
include("testNLsolve.jl")
include("testCommonConvWrapper.jl")

include("testApproxConv.jl")
include("testBasicForwardConvolve.jl")
include("testUseMsgLikelihoods.jl")
include("testDefaultDeconv.jl")

include("testPartialFactors.jl")
include("testPartialPrior.jl")
include("testpartialconstraint.jl")
include("testPartialNH.jl")
include("testMixturePrior.jl")

include("testStateMachine.jl")
include("testBasicCSM.jl")
include("testCliqueFactors.jl")
include("testCcolamdOrdering.jl")
include("testBasicGraphs.jl")
include("testJointEnforcement.jl")
include("testHasPriors913.jl")
include("testInitVariableOrder.jl")
include("testTreeMessageUtils.jl")
include("testCSMMonitor.jl")
include("testExpXstroke.jl")
include("testBasicRecycling.jl")
include("testSkipUpDown.jl")
include("testlocalconstraintexamples.jl")
include("testManualInit.jl")
include("testBasicTreeInit.jl")
include("testSolveOrphanedFG.jl")
include("testSolveSetPPE.jl")
include("testSolveKey.jl")
end

if TEST_GROUP in ["all", "test_cases_group"]
include("testnullhypothesis.jl") 
include("testVariousNSolveSize.jl")
include("testExplicitMultihypo.jl")
include("TestCSMMultihypo.jl")
include("testCalcFactorHypos.jl")
include("testMultimodal1D.jl")
include("testMultihypoAndChain.jl")
include("testMultithreaded.jl")
include("testmultihypothesisapi.jl")
include("fourdoortest.jl")
include("testCircular.jl")
include("testMixtureLinearConditional.jl")
include("testFluxModelsDistribution.jl")
include("testAnalysisTools.jl")
include("testDERelative.jl")

include("testBasicParametric.jl")
include("testMixtureParametric.jl")

# dont run test on ARM, as per issue #527
if Base.Sys.ARCH in [:x86_64;]
  include("testTexTreeIllustration.jl")
end

include("testMultiprocess.jl")
include("testDeadReckoningTether.jl")
end


#
