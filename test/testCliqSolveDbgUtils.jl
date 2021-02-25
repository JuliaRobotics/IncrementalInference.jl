using IncrementalInference
using Test

@testset "Test SolveCliqueUp!" begin
##
N=8
fg = generateCanonicalFG_lineStep(N; 
                                  graphinit=false,
                                  poseEvery=1, 
                                  landmarkEvery=N+1, 
                                  posePriorsAt=[0],
                                  landmarkPriorsAt=[], 
                                  sightDistance=N+1)

deleteFactor!.(fg, [Symbol("x$(i)lm0f1") for i=1:(N-1)])

ensureAllInitialized!(fg)

tree = buildTreeReset!(fg)

a,b = solveCliqUp!(fg, tree, 2)
##
end
