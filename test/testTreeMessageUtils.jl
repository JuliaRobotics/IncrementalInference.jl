
using Test
using IncrementalInference

@testset "Testing tree message utils" begin

N=8
fg = generateGraph_LineStep(N; 
                            graphinit=false,
                            poseEvery=1, 
                            landmarkEvery=N+1, 
                            posePriorsAt=[0],
                            landmarkPriorsAt=[], 
                            sightDistance=N+1)

deleteFactor!.(fg, [Symbol("x$(i)lm0f1") for i=1:(N-1)])

# fixed eliminationOrder for repeatability
eliminationOrder = [:x3, :x8, :x5, :x1, :x6, :lm0, :x7, :x4, :x2, :x0]
smtasks = Task[]
tree = solveTree!(fg; smtasks=smtasks, eliminationOrder=eliminationOrder)

allmsgs = getTreeCliqUpMsgsAll(tree)
#TODO better test but for now spot check a few keys
@test issetequal(collect(2:8), keys(allmsgs)) 
belief2 = allmsgs[2].belief
@test issetequal([:x0, :x4], keys(belief2)) 


stackedmsgs = stackCliqUpMsgsByVariable(tree, allmsgs)
#everything except leave frontals
allvars = [:lm0, :x0, :x2, :x4, :x6, :x7]
@test issetequal(allvars, keys(stackedmsgs)) 
@test length(stackedmsgs[:x0]) == 3
x4stack =  stackedmsgs[:x4]
cliq_depths = map(v->v.cliqId[]=>v.depth, x4stack)
@test issetequal(cliq_depths, [4 => 2, 6 => 3, 2 => 1, 8 => 1])

c3 = getClique(tree, IIF.CliqueId(3))
c7 = getClique(tree, IIF.CliqueId(7))

@test IIF.compare(IIF.getMessageUpRx(c3)[7], IIF.getMessageUpTx(c7))
@test IIF.compare(IIF.getMessageDownTx(c3), IIF.getMessageDownRx(c7))

end