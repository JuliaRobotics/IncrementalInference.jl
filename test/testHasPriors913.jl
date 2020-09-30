using IncrementalInference
using Test


@testset "Only do UPWARD_COMMON message likelihoods if hasPriors" begin


fg = generateCanonicalFG_lineStep(4; 
                                  poseEvery=1, 
                                  landmarkEvery=5, 
                                  posePriorsAt=[],
                                  landmarkPriorsAt=[], 
                                  sightDistance=5,
                                  solverParams=SolverParams(algorithms=[:default, :parametric]))

deleteFactor!.(fg, [Symbol("x$(i)lm0f1") for i=1:3])

#force wrong init
prpo = Prior(Normal(5, 0.01))
addFactor!(fg, [:x0], prpo)
ensureAllInitialized!(fg)
deleteFactor!(fg, :x0f1)

# now the correct prior
prpo = Prior(Normal(0, 0.01))
addFactor!(fg, [:x0], prpo)

fg.solverParams.useMsgLikelihoods = true

smtasks = Task[]
tree, smt, hists = solveTree!(fg; smtasks=smtasks, verbose=true, timeout=30);
# tree, smt, hists = solveTree!(fg; smtasks, verbose=true, timeout=20, recordcliqs=ls(fg));


for i = 0:4
  # btd = getBelief()
  @show i getPPE(getVariable(fg, Symbol("x$i"))).suggested
  @error("Must enable hasPriors test")
    # @test isapprox(mean(getKDEfit(btd,distribution=Normal)), i; atol=0.1)
end 


end





#
