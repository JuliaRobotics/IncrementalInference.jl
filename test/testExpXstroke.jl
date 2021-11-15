
##

using IncrementalInference
using Test

##

@testset "test endless cycle case, issue #754" begin

##

# testgraph from issue #754
fg = generateGraph_LineStep(5; 
                            poseEvery=1, 
                            landmarkEvery=5, 
                            posePriorsAt=[0,2], 
                            sightDistance=4)
#
getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true
getSolverParams(fg).useMsgLikelihoods = true


##

# tree = buildTreeReset!(fg)
# drawTree(tree, show=true)


##

# ENV["JULIA_DEBUG"] = :csm_4

smtasks = Task[]
hist = IIF.solveTree!(fg; smtasks=smtasks); #, recordcliqs=ls(fg));


##

# hists = fetchCliqHistoryAll!(smtasks);
# printCSMHistorySequential(hists)
# printCSMHistoryLogical(hists)

# ##

# fnc_ = repeatCSMStep!(hists, 4, 6)


##


for var in sortDFG(ls(fg))
  sppe = getVariable(fg,var) |> getPPE |> IIF.getPPESuggested
  println("Testing ", var,": ", sppe)
  @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.2)
end


##


# linear octo 
N=8
fg = generateGraph_LineStep(N; 
                            graphinit=false,
                            poseEvery=1, 
                            landmarkEvery=N+1, 
                            posePriorsAt=[0],
                            landmarkPriorsAt=[], 
                            sightDistance=N+1)

deleteFactor!.(fg, [Symbol("x$(i)lm0f1") for i=1:(N-1)])

getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true
getSolverParams(fg).useMsgLikelihoods = true

smtasks = Task[]
tree = IIF.solveTree!(fg; smtasks=smtasks);


##

# hists = fetchCliqHistoryAll!(smtasks);
# printCSMHistorySequential(hists, 2=>1:50)
# printCSMHistoryLogical(hists)

# ##

# csmc_ = repeatCSMStep!(hists, 2, 14)


##

for var in sortDFG(ls(fg))
    sppe = getVariable(fg,var) |> getPPE |> IIF.getPPESuggested
    println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.4)
end


##

# Larger graph
fg = generateGraph_LineStep(15; 
                            poseEvery=1, 
                            landmarkEvery=3, 
                            posePriorsAt=[0,7,12],
                            landmarkPriorsAt=[0,3], 
                            sightDistance=2)

getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true
getSolverParams(fg).useMsgLikelihoods = true

smtasks = Task[]
tree = IIF.solveTree!(fg; smtasks=smtasks);

for var in sortDFG(ls(fg))
    sppe = getVariable(fg,var) |> getPPE |> IIF.getPPESuggested
    println("Testing ", var,": ", sppe)
    s = findfirst(r"\d", string(var))[1]
    @test isapprox(sppe[1], parse(Int,string(var)[s:end]), atol=0.45)
end

##

end