using IncrementalInference
using Test

## a new factor that is broken
struct BrokenFactor{T<: SamplableBelief} <: AbstractManifoldMinimize
    Z::T
end

IIF.getManifold(::BrokenFactor) = TranslationGroup(1)
function IIF.getSample(cf::CalcFactor{<:BrokenFactor})
  return rand(cf.factor.Z, 1)
end

function (s::CalcFactor{<:BrokenFactor})(z,
                                        wxi,
                                        wxj)
    #
    error("User factor has a bug.")
end

# # FIXME consolidate with CalcFactor according to #467
# function (s::BrokenFactor{<:IIF.ParametricTypes})(X1::AbstractArray{<:Real},
#                                                   X2::AbstractArray{<:Real};
#                                                   userdata::Union{Nothing,FactorMetadata}=nothing )
#     error("User factor has a bug -- USE NEW CalcFactor API INSTEAD, v0.21.")
# end

##

@testset "Test CSM monitor/watchdog on errors" begin

##

# create a factor graph
fg = generateGraph_LineStep(10; 
                            poseEvery=1, 
                            landmarkEvery=5, 
                            posePriorsAt=[0], 
                            sightDistance=4,
                            solverParams=SolverParams(algorithms=[:default, :parametric]))
#

initAll!(fg)

##

#TODO test FSM watchdog
# add a broken factor - mid
addFactor!(fg, [:x9, :lm10], BrokenFactor(Normal()); graphinit=false)
smtasks = Task[]
@test_throws CompositeException tree = IIF.solveTree!(fg; smtasks=smtasks);
sleep(0.1)

## Test parametric solve also



addFactor!(fg, [:x9, :lm10], BrokenFactor(Normal()); graphinit=false)
##
# IIF.solveTree!(fg; smtasks=smtasks, algorithm = :parametric)
##
@test_throws CompositeException tree2 = IIF.solveTree!(fg; smtasks=smtasks, algorithm = :parametric)
sleep(0.1)

deleteFactor!(fg, :x9lm10f2)

## add a broken factor - leave

addFactor!(fg, [:x10, :lm10], BrokenFactor(Normal()); graphinit=false)
@test_throws CompositeException tree2 = IIF.solveTree!(fg; smtasks=smtasks, algorithm = :parametric)
sleep(0.1)

deleteFactor!(fg, :x10lm10f2)

## add a broken factor - root
addFactor!(fg, [:x7, :lm10], BrokenFactor(Normal()); graphinit=false)
@test_throws CompositeException tree2 = IIF.solveTree!(fg; smtasks=smtasks, algorithm = :parametric)


##

end


@testset "test CSM debug options work" begin

## create a factor graph

fg = generateGraph_LineStep(10; 
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
csmc_ = repeatCSMStep!(hists, 1, 1);


##

end