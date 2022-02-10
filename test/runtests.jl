using Test

# temporarily moved to start (for debugging)
include("testDefaultDeconv.jl")
include("testSpecialEuclidean2Mani.jl")

include("testFactorGradients.jl")
include("testpartialconstraint.jl")
include("testPartialNH.jl")

include("testSphereMani.jl")
include("testSpecialOrthogonalMani.jl")

include("testDistributionsGeneric.jl")
include("testHeatmapGridDensity.jl")
include("testCliqSolveDbgUtils.jl")
include("TestModuleFunctions.jl")
include("testApproxConv.jl")
include("testCompareVariablesFactors.jl")
include("typeReturnMemRef.jl")
include("basicGraphsOperations.jl")
include("testMixturePrior.jl")
include("testGradientUtils.jl")
include("testPartialFactors.jl")
include("testPartialPrior.jl")
include("testSaveLoadDFG.jl")
include("testJunctionTreeConstruction.jl")
include("testBayesTreeiSAM2Example.jl")
include("testTreeFunctions.jl")

#FIXME fails on MetaBayesTree
include("testTreeSaveLoad.jl")

include("testSpecialSampler.jl")
include("saveconvertertypes.jl")
include("testgraphpackingconverters.jl")
include("testNLsolve.jl")
include("testCommonConvWrapper.jl")
include("testBasicForwardConvolve.jl")
include("testFactorMetadata.jl")
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
include("testBasicTreeInit.jl")
include("testSolveOrphanedFG.jl")
include("testSolveSetPPE.jl")
include("testSolveKey.jl")
include("testEuclidDistance.jl")
include("priorusetest.jl")
include("testnullhypothesis.jl") 
include("testVariousNSolveSize.jl")
include("testExplicitMultihypo.jl")
include("TestCSMMultihypo.jl")
include("testMultihypoFMD.jl")
include("testMultiHypo2Door.jl")
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



#
