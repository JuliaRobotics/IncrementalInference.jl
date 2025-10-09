using IncrementalInference
using Test


##

@testset "Only do UPWARD_COMMON message likelihoods if hasPriors" begin

##


fg = generateGraph_LineStep(4; 
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
initAll!(fg)
deleteFactor!(fg, :x0f1)

# now the correct prior
prpo = Prior(Normal(0, 0.01))
addFactor!(fg, [:x0], prpo)

fg.solverParams.useMsgLikelihoods = true

smtasks = Task[]
tree = solveTree!(fg; smtasks=smtasks, verbose=true, timeout=30);

@warn("hasPriors test needs multiple solves")
tree = solveTree!(fg);
tree = solveTree!(fg);
# tree = solveTree!(fg; smtasks, verbose=true, timeout=20, recordcliqs=ls(fg));

for i = 0:4
  ppe = calcMeanMaxSuggested(fg, Symbol("x$i"), :default).suggested[1]
  @show i ppe
  @test isapprox(ppe, i; atol=0.7)
end 


##

end





#
