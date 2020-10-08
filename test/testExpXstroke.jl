using IncrementalInference
using Test


@testset "test endless cycle case, issue #754" begin

# testgraph from issue #754
fg = generateCanonicalFG_lineStep(5; 
                                  poseEvery=1, 
                                  landmarkEvery=5, 
                                  posePriorsAt=[0,2], 
                                  sightDistance=4)
                                  
getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true

smtasks = Task[]
tree, smt, hist = IIF.solveTree_X!(fg; smtasks=smtasks);

for var in sortDFG(ls(fg))
  sppe = getVariable(fg,var) |> getPPE |> IIF.getSuggestedPPE
  println("Testing ", var,": ", sppe)
  @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.1)
end


# linear octo 
N=8
fg = generateCanonicalFG_lineStep(N; 
                                  graphinit=false,
                                  poseEvery=1, 
                                  landmarkEvery=N+1, 
                                  posePriorsAt=[0],
                                  landmarkPriorsAt=[], 
                                  sightDistance=N+1)

deleteFactor!.(fg, [Symbol("x$(i)lm0f1") for i=1:(N-1)])

getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true

smtasks = Task[]
tree, smt, hists = IIF.solveTree_X!(fg; smtasks=smtasks);

for var in sortDFG(ls(fg))
    sppe = getVariable(fg,var) |> getPPE |> IIF.getSuggestedPPE
    println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.15)
end


# Larger graph
fg = generateCanonicalFG_lineStep(15; 
                                  poseEvery=1, 
                                  landmarkEvery=3, 
                                  posePriorsAt=[0,7,12],
                                  landmarkPriorsAt=[0,3], 
                                  sightDistance=2)

getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true

smtasks = Task[]
tree, smt, hists = IIF.solveTree_X!(fg; smtasks=smtasks);

for var in sortDFG(ls(fg))
    sppe = getVariable(fg,var) |> getPPE |> IIF.getSuggestedPPE
    println("Testing ", var,": ", sppe)
    s = findfirst(r"\d", string(var))[1]
    @test isapprox(sppe[1], parse(Int,string(var)[s:end]), atol=0.2)
end


end