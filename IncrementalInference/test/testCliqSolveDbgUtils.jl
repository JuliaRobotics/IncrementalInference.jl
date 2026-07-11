using IncrementalInference
using Test

##

@testset "Test solveCliqueUp! and solveCliqDown!" begin

##

N=8
fg = generateGraph_LineStep(N; 
                            graphinit=false,
                            poseEvery=1, 
                            landmarkEvery=N+1, 
                            posePriorsAt=[0],
                            landmarkPriorsAt=[], 
                            sightDistance=N+1)
#

deleteFactor!.(fg, [Symbol("x$(i)lm0f1") for i=1:(N-1)])

# test the initAll! separately anyway
initAll!(fg)

tree = buildTreeReset!(fg)

# for debuggin use 
# ENV["JULIA_DEBUG"] = :csm_2

# solve clique up tests
a,b = solveCliqUp!(fg, tree, 2)
a,b = solveCliqUp!(fg, tree, 2; recordcliq = true) 
@test length(a) > 0

# solve clique down tests 
a,b = solveCliqDown!(fg, tree, 2)
a,b = solveCliqDown!(fg, tree, 2; recordcliq = true) 
@test length(a) > 0

##

end



