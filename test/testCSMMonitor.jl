using IncrementalInference
using Test

## a new factor that is broken
struct BrokenFactor{T<: SamplableBelief} <: AbstractRelativeRoots
    Z::T
end

IncrementalInference.getSample(cf::CalcFactor{<:BrokenFactor}, N::Int=1) = (reshape(rand(cf.factor.Z, N),:,N), )

function (s::CalcFactor{<:BrokenFactor})(z,
                                         wxi,
                                         wxj)
    #
    error("User factor has a bug.")
end

# FIXME consolidate with CalcFactor according to #467
function (s::BrokenFactor{<:IIF.ParametricTypes})(X1::AbstractArray{<:Real},
                                                  X2::AbstractArray{<:Real};
                                                  userdata::Union{Nothing,FactorMetadata}=nothing )
    error("User factor has a bug.")
end

@testset "Test CSM monitor/watchdog on errors" begin

# create a factor graph
fg = generateCanonicalFG_lineStep(10; 
                                  poseEvery=1, 
                                  landmarkEvery=5, 
                                  posePriorsAt=[0], 
                                  sightDistance=4,
                                  solverParams=SolverParams(algorithms=[:default, :parametric]))
                            
ensureAllInitialized!(fg)

#TODO test FSM watchdog
# add a broken factor - mid
addFactor!(fg, [:x9, :lm10], BrokenFactor(Normal()); graphinit=false)
smtasks = Task[]
@test_throws CompositeException tree, smt, hist = IIF.solveTree!(fg; smtasks=smtasks);
sleep(0.1)

# Test parametric solve also
addFactor!(fg, [:x9, :lm10], BrokenFactor(Normal()); graphinit=false)
@test_throws CompositeException tree2, smt, hist = IIF.solveTree!(fg; smtasks=smtasks, algorithm = :parametric)
sleep(0.1)

deleteFactor!(fg, :x9lm10f2)
# add a broken factor - leave
addFactor!(fg, [:x10, :lm10], BrokenFactor(Normal()); graphinit=false)
@test_throws CompositeException tree2, smt, hist = IIF.solveTree!(fg; smtasks=smtasks, algorithm = :parametric)
sleep(0.1)

deleteFactor!(fg, :x10lm10f2)
# add a broken factor - root
addFactor!(fg, [:x7, :lm10], BrokenFactor(Normal()); graphinit=false)
@test_throws CompositeException tree2, smt, hist = IIF.solveTree!(fg; smtasks=smtasks, algorithm = :parametric)

end




@testset "test CSM debug options work" begin


## create a factor graph

fg = generateCanonicalFG_lineStep(10; 
                                  poseEvery=1, 
                                  landmarkEvery=5, 
                                  posePriorsAt=[0], 
                                  sightDistance=4 )
#                            

smtasks = Task[];
solveTree!(fg, smtasks=smtasks, recordcliqs=ls(fg));


## make sure we fetch the results


hists = fetchCliqHistoryAll!(smtasks);

printCSMHistoryLogical(hists)

printCSMHistorySequential(hists)

## test CSM resolve steps

@info "test repeatCSMStep"
csmc_ = repeatCSMStep!(hists, 1, 1)


##

end