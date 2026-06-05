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
    @testset "Temporary Debug Group (most development activity, fail fast)" begin
            
      include("testSpecialSampler.jl") # sometimes BoundsError, suspect need resample to N step, EvaluFactor.jl:179
        include("testDERelative.jl") # FIX BoundsError Ln332 cf._legacyParams[k][i], [100] of 1..99
      
      include("testEuclidDistance.jl") # test_broken
      include("testMultiHypo3Door.jl") # FIX, slow, weak numerics
      include("testExpXstroke.jl") # FIX, init bw issue, addLikelihoodsDifferentialCHILD! LN327
      include("testCircular.jl") # FIX
      include("testFluxModelsDistribution.jl") # FIX
      include("testMixturePrior.jl") # FIX, serde structutils issue with BinarTruckFixedDepth?
      include("testMixtureLinearConditional.jl") # FIX, use HomotopyDensity as replacement for Mixture
      include("testMixtureParametric.jl") #FIXME parametric mixtures #1787
      @test_broken error("testSphereMani.jl broken") # include("testSphereMani.jl") # FIXME

      include("testSpecialOrthogonalMani.jl") # FIX, stateLabel :parametric not found

      include("testSpecialEuclidean2Mani.jl") # FIX, parallel_transport_curvature_2nd_lie not defined for this case
      include("testpartialconstraint.jl") # FIX, big numerical fail
      include("testPartialNH.jl") # FIX
      include("testCcolamdOrdering.jl") # FIX

      include("testMultimodal1D.jl") # FIX numerics, skipped a test
      # include("testMultiprocess.jl")

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
    @testset "Basic Functional Group (stable functional tests)" begin

      include("testBasicParametric.jl") # FIX, access undef ref


      # start as basic as possible and build from there
      include("typeReturnMemRef.jl")
      include("testDistributionsGeneric.jl")
      include("basicGraphsOperations.jl")

      #FIXME fails on MetaBayesTree
      include("testTreeSaveLoad.jl")

      # test convolution functions
      include("testApproxConv.jl")
      include("testBasicForwardConvolve.jl")
      
      include("testCliqSolveDbgUtils.jl")

      include("testBasicManifolds.jl")

      include("testJunctionTreeConstruction.jl")
      include("testBayesTreeiSAM2Example.jl")


      include("testCommonConvWrapper.jl") # skipped test with just ::Float point type 

      include("testStateMachine.jl")
      include("testBasicCSM.jl")
      include("testCliqueFactors.jl")

      include("testCliqueTreesOrderings.jl")
      include("testlocalconstraintexamples.jl")

      # dont run test on ARM, as per issue #527
      if Base.Sys.ARCH in [:x86_64;]
        include("testTexTreeIllustration.jl")
      end

      include("testManualInit.jl")
      include("testSolveOrphanedFG.jl")
      include("testSolveKey.jl")

      include("testJointEnforcement.jl")
      include("testPartialFactors.jl")
      include("testCalcFactorHypos.jl")

      include("testnullhypothesis.jl")
      include("testExplicitMultihypo.jl")
      include("TestCSMMultihypo.jl")

      include("testMultithreaded.jl")
      include("testmultihypothesisapi.jl")
      include("testAnalysisTools.jl")
      include("testVariousNSolveSize.jl")

      include("testDeadReckoningTether.jl")

      include("testHeatmapGridDensity.jl")

      # refac AMP v0.15, these functionals are medium slow
      include("testPartialPrior.jl") # slowish
      include("testDefaultDeconv.jl") # slowish
      include("testUseMsgLikelihoods.jl") # slowish
      include("testTreeFunctions.jl") # slower
      include("testInitVariableOrder.jl") # slowish
      include("testTreeMessageUtils.jl") # slowish
      include("testMultihypoAndChain.jl") # slower, weak numerics

    end
  end # basic_functional_group


  if TEST_GROUP in ["all", "test_cases_group"]
    @testset "Test Cases Group (offload concurrent CI, slow running jobs)" begin
      
      # refac AMP v0.15, these tests are very slow
      include("fourdoortest.jl") # slowish
      include("testBasicGraphs.jl") # slow, weak numerics
      include("priorusetest.jl") # slow, many skips
      include("testHasPriors913.jl") # slow

      include("testBasicTreeInit.jl") # slow
      include("testBasicRecycling.jl") # slow
      include("testSkipUpDown.jl") # slow
      include("testCSMMonitor.jl") # slow

    end
  end # test_cases_group

end # all IIF tests

#
