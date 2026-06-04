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
    @testset "Temporary Debug Group (where most development work is happening)" begin
      include("testBasicGraphs.jl") # NUMERICAL
      
      include("testEuclidDistance.jl") # test_broken
      include("testSpecialSampler.jl") # FIX, bounds error
      include("testDERelative.jl") # FIX BoundsError Ln332 cf._legacyParams[k][i], [100] of 1..99
      include("testHeatmapGridDensity.jl") # numerical test broken
      include("testMultiHypo3Door.jl") # FIX numerical
      include("testExpXstroke.jl") # FIX, init bw issue, addLikelihoodsDifferentialCHILD! LN327
      include("testCircular.jl") # FIX
      include("testMultihypoAndChain.jl") # numerical issues
      include("testFluxModelsDistribution.jl")
      include("testMixturePrior.jl") # FIX, serde structutils issue with BinarTruckFixedDepth?
      include("testMixtureLinearConditional.jl") # FIX, use HomotopyDensity as replacement for Mixture
      include("testMixtureParametric.jl") #FIXME parametric mixtures #1787
      @test_broken error("testSphereMani.jl broken") # include("testSphereMani.jl") # FIXME

      # gradient / jacobian tests
      #include("manifolds/manifolddiff.jl")
      #include("manifolds/factordiff.jl")
      @error "Gradient tests must be updated and restored for new ccw.varValsAll[]"
      #include("testGradientUtils.jl")
      #include("testFactorGradients.jl")

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
  end # tmp_debug_group


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

      include("priorusetest.jl") # slow, many skips

      include("testCliqSolveDbgUtils.jl")
      include("testUseMsgLikelihoods.jl")

      include("testBasicManifolds.jl")

      include("testJunctionTreeConstruction.jl")
      include("testBayesTreeiSAM2Example.jl")
      include("testTreeFunctions.jl")

      include("testCommonConvWrapper.jl") # skipped test with just ::Float point type 

      include("testStateMachine.jl")
      include("testBasicCSM.jl") # ??, NLLsolver options has no field .defaultNumKernels
      include("testCliqueFactors.jl")
      include("testCcolamdOrdering.jl")
      include("testCliqueTreesOrderings.jl")
      include("testHasPriors913.jl")
      include("testInitVariableOrder.jl")
      include("testTreeMessageUtils.jl")
      include("testBasicRecycling.jl") # slow
      include("testSkipUpDown.jl") # slow
      include("testlocalconstraintexamples.jl")

      include("testManualInit.jl")
      include("testBasicTreeInit.jl")
      include("testSolveOrphanedFG.jl")
      include("testSolveKey.jl")

      include("testJointEnforcement.jl")
      include("testPartialFactors.jl")
      include("testCalcFactorHypos.jl")

      include("testnullhypothesis.jl")
      include("testExplicitMultihypo.jl")
      include("TestCSMMultihypo.jl")

      include("testPartialPrior.jl")
      include("testMultithreaded.jl")
      include("testmultihypothesisapi.jl")
      include("testAnalysisTools.jl")
      include("testVariousNSolveSize.jl")
      include("fourdoortest.jl")

      include("testDeadReckoningTether.jl") # FIX convert vector to set, VariableDFG

      # dont run test on ARM, as per issue #527
      if Base.Sys.ARCH in [:x86_64;]
        include("testTexTreeIllustration.jl")
      end
    end
  end # basic_functional_group


  if TEST_GROUP in ["all", "test_cases_group"]
    @testset "Test Cases Group (offload concurrent CI of slow running jobs)" begin
      ## WORK IN PROGRESS

      include("testCSMMonitor.jl") # slow, FieldError: type IncrementalInference.NLLSSolver has no field `defaultNumKernels`; IncrementalInference.NLLSSolver has no fields at all. GraphInit.jlLn14 prepareState!
      include("testSpecialOrthogonalMani.jl") # FIX

      include("testSpecialEuclidean2Mani.jl") # TBD
      include("testpartialconstraint.jl") # FIX
      include("testPartialNH.jl") # FIX
      include("testBasicParametric.jl") # FIX

      include("testMultimodal1D.jl") # FIX numerics, skipped a test
      # include("testMultiprocess.jl")
    end
  end # test_cases_group

end # all IIF tests

#
