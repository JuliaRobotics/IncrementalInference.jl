# addprocs(3)
using Test
# using Compat
# using IncrementalInference


include("testBasicParametric.jl")

include("TestModuleFunctions.jl")

include("testStateMachine.jl")

include("testCompareVariablesFactors.jl")

include("typeReturnMemRef.jl")

include("basicGraphsOperations.jl")

include("testMixturePrior.jl")

include("testPartialFactors.jl")

include("testBayesTreeiSAM2Example.jl")

include("testSpecialSampler.jl")

include("testSaveLoadDFG.jl")

include("testJunctionTreeConstruction.jl")

#FIXME fails on MetaBayesTree
include("testTreeSaveLoad.jl")

include("saveconvertertypes.jl")

include("testgraphpackingconverters.jl")

include("testNLsolve.jl")

include("testNumericRootGenericRandomized.jl")

include("testCommonConvWrapper.jl")

include("testBasicForwardConvolve.jl")

include("testFactorMetadata.jl")

include("testBasicCSM.jl")

include("testCliqueFactors.jl")

include("testCcolamdOrdering.jl")

include("testBasicGraphs.jl")

# TODO old names should be removed, like Odo, Obsv2
include("testlocalconstraintexamples.jl")

include("testBasicTreeInit.jl")

include("testSolveOrphanedFG.jl")

include("testSolveSetPPE.jl")

# include("priorusetest.jl")

include("testVariousNSolveSize.jl")

include("testExplicitMultihypo.jl")

include("TestCSMMultihypo.jl")

include("testMultiHypo2Door.jl")

include("testMultimodal1D.jl")

include("testMultithreaded.jl")

include("testpartialconstraint.jl")

include("testnullhypothesis.jl")

include("testmultihypothesisapi.jl")

include("fourdoortest.jl")

include("testAnalysisTools.jl")

# dont run test on ARM, as per issue #527
if Base.Sys.ARCH in [:x86_64;]
  include("testTexTreeIllustration.jl")
end










#
