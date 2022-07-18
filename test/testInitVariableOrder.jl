using Test
using IncrementalInference

@testset "Testing getCliqVarInitOrderUp" begin
  
fg = generateGraph_LineStep(3; 
                            poseEvery=1, 
                            landmarkEvery=5, 
                            posePriorsAt=[0],
                            landmarkPriorsAt=[0,2], 
                            sightDistance=3)                                  

fg.solverParams.useMsgLikelihoods = true
# addVariable!(subfg, :x0, Con)

cliqInitOrder = getCliqVarInitOrderUp(fg) 
#TODO variable order can change, but should be more stable using OrderedDict. Testing for both.
# maybe remove cliqInitOrder == [:x0, :lm0, :x3, :x2, :x1] test in v0.31
@test cliqInitOrder == [:x0, :lm0, :x3, :x2, :x1] || cliqInitOrder == [:x0, :lm0, :x3, :x1, :x2]

solveTree!(fg)

# construct a message on :x3
msg = LikelihoodMessage(status=IIF.UPSOLVED)
seps = [:x3]
for vid in seps
  var = DFG.getVariable(fg, vid)
  if isInitialized(var)
    msg.belief[var.label] = TreeBelief(var)
  end
end
addMsgFactors!(fg, msg, IIF.UpwardPass)

# less sensitive to order given msgjoint message updates
@test intersect(getCliqVarInitOrderUp(fg), [:x3, :x0, :lm0, :x2, :x1]) |> length == 5

# construct a message on :lm0,:x2
msg = LikelihoodMessage(status=IIF.UPSOLVED, hasPriors=false)
seps = [:lm0, :x2]
for vid in seps
    var = DFG.getVariable(fg, vid)
    if isInitialized(var)
        msg.belief[var.label] = TreeBelief(var)
    end
end
addMsgFactors!(fg, msg, IIF.UpwardPass)

@warn("disabled getCliqVarInitOrderUp test until #1010 is closed")
# @test getCliqVarInitOrderUp(fg) == [:x3, :x0, :lm0, :x1, :x2]


# test order with mixture prior #998
cv = 3.0
doorPrior = Mixture(Prior,
                    [Normal(-100,cv);Normal(0,cv);Normal(100,cv);Normal(300,cv)],
                    [1/4;1/4;1/4;1/4] )

                            
fg = initfg()
fg.solverParams.graphinit = false
addVariable!(fg, :x0, ContinuousScalar)

addFactor!(fg, [:x0], doorPrior)

addVariable!(fg, :x1, ContinuousScalar)

addFactor!(fg,[:x0; :x1], LinearRelative( Normal(0.0,1.0)))

@test lsfPriors(fg) == [:x0f1]

@test getCliqVarInitOrderUp(fg) == [:x0, :x1]

end