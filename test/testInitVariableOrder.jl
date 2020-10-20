using Test
using IncrementalInference

@testset "Testing getCliqVarInitOrderUp" begin
  
fg = generateCanonicalFG_lineStep(3; 
                                  poseEvery=1, 
                                  landmarkEvery=5, 
                                  posePriorsAt=[0],
                                  landmarkPriorsAt=[0,2], 
                                  sightDistance=3)                                  

fg.solverParams.useMsgLikelihoods = true
# addVariable!(subfg, :x0, Con)

@test getCliqVarInitOrderUp(fg) == [:x0, :lm0, :x3, :x2, :x1]

solveTree!(fg)

# construct a message on :x3
msg = LikelihoodMessage(status=:UPSOLVED)
seps = [:x3]
for vid in seps
  var = DFG.getVariable(fg, vid)
  if isInitialized(var)
    msg.belief[var.label] = TreeBelief(var)
  end
end
addMsgFactors!(fg, msg, IIF.UpwardPass)

@test getCliqVarInitOrderUp(fg) == [:x3, :x0, :lm0, :x2, :x1]

# construct a message on :lm0,:x2
msg = LikelihoodMessage(status=:UPSOLVED, hasPriors=false)
seps = [:lm0, :x2]
for vid in seps
    var = DFG.getVariable(fg, vid)
    if isInitialized(var)
        msg.belief[var.label] = TreeBelief(var)
    end
end
addMsgFactors!(fg, msg, IIF.UpwardPass)

@test getCliqVarInitOrderUp(fg) == [:x3, :x0, :lm0, :x1, :x2]

end