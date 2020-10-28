
using Test
using IncrementalInference

@testset "Skip downsolve" begin
N=6
fg = generateCanonicalFG_lineStep(N; 
                                  graphinit=false,
                                  poseEvery=1, 
                                  landmarkEvery=N+1, 
                                  posePriorsAt=[0],
                                  landmarkPriorsAt=[], 
                                  sightDistance=N+1)

deleteFactor!.(fg, [Symbol("x$(i)lm0f1") for i=1:(N-1)])

getSolverParams(fg).useMsgLikelihoods = true
getSolverParams(fg).downsolve = false

smtasks = Task[]
tree, smt, hists = solveTree!(fg; smtasks=smtasks, recordcliqs=[:x4]);

# See if downsolve was called
@test !(IIF.solveDown_StateMachine in getindex.(hists[2], 3))

#test if values are still correct
for var in sortDFG(ls(fg))
    sppe = getVariable(fg,var) |> getPPE |> IIF.getPPESuggested
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.2)
end


#now downsolve only
getSolverParams(fg).upsolve = false
getSolverParams(fg).downsolve = true

tree, smt, hists = solveTree!(fg; smtasks=smtasks, recordcliqs=[:x4]);
# See if upsolved was called
@test !(IIF.solveUp_StateMachine in getindex.(hists[2], 3))

#test if values are still correct
for var in sortDFG(ls(fg))
    sppe = getVariable(fg,var) |> getPPE |> IIF.getPPESuggested
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.2)
end

# recycled downsolve only
tree, smt, hists = solveTree!(fg, tree; smtasks=smtasks, recordcliqs=[:x4]);
@test !(IIF.solveUp_StateMachine in getindex.(hists[2], 3))



end