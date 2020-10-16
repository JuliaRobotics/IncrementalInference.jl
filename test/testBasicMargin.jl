using Test
using IncrementalInference

@testset "Testing basic marginalization" begin
# linear hex 
N=6
fg = generateCanonicalFG_lineStep(N; 
                                  graphinit=false,
                                  poseEvery=1, 
                                  landmarkEvery=N+1, 
                                  posePriorsAt=[0],
                                  landmarkPriorsAt=[], 
                                  sightDistance=N+1)

deleteFactor!.(fg, [Symbol("x$(i)lm0f1") for i=1:(N-1)])

# tree = resetBuildTree!(fg, show=true, drawpdf=true)
tree, smtasks, hists = solveTree!(fg)

for var in sortDFG(ls(fg))
    sppe = getVariable(fg,var) |> getPPE |> IIF.getSuggestedPPE
    println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.1)
end


defaultFixedLagOnTree!(fg, 6)
# Back up data from these two poses so we can compare them once we solve again.
lm0 = deepcopy(getVal(fg, :lm0))
X0 = deepcopy(getVal(fg, :x0))
X1 = deepcopy(getVal(fg, :x1))

fifoFreeze!(fg)

tree, smtasks, hists = solveTree!(fg)
for var in sortDFG(ls(fg))
    sppe = getVariable(fg,var) |> getPPE |> IIF.getSuggestedPPE
    println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.15)
end


lm0cmp = deepcopy(getVal(fg, :lm0))
X0cmp = deepcopy(getVal(fg, :x0))
X1cmp = deepcopy(getVal(fg, :x1))

@test X0 == X0cmp #Frozen
@test lm0 == lm0cmp #Frozen
@test X1 != X1cmp #Recalculated

deleteVariable!(fg, :x0)

addVariable!.(fg, [Symbol("x$i") for i = 7:9], ContinuousScalar)
addFactor!(fg, [:x6,:x7], LinearRelative(Normal(1.0, 0.1)))
addFactor!(fg, [:x7,:x8], LinearRelative(Normal(1.0, 0.1)))
addFactor!(fg, [:x8,:x9], LinearRelative(Normal(1.0, 0.1)))
addFactor!(fg, [:lm0, :x9], LinearRelative(Normal(9,0.1)))

tree, smtasks, hists = solveTree!(fg)

for var in sortDFG(ls(fg))
    sppe = getVariable(fg,var) |> getPPE |> IIF.getSuggestedPPE
    println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.15)
end

@test lm0 == getVal(fg, :lm0) #Still Frozen
@test X1cmp == getVal(fg, :x1) #also now Frozen

unfreezeVariablesAll!(fg)

# freeze again
defaultFixedLagOnTree!(fg, 9)

tree, smtasks, hists = solveTree!(fg)

@test lm0 == getVal(fg, :lm0) #Still Frozen
@test X1cmp != getVal(fg, :x1) #not frozen

# freeze 2,4,6 to all marginalize clique 2
setfreeze!(fg, [:x2, :x4, :x6])
tree, smtasks, hists = solveTree!(fg)

@test areCliqVariablesAllMarginalized(fg, tree.cliques[2])

tree, smtasks, hists = solveTree!(fg, tree)#, recordcliqs=ls(fg));
for var in sortDFG(ls(fg))
    sppe = getVariable(fg,var) |> getPPE |> IIF.getSuggestedPPE
    println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.15)
end

X1 = deepcopy(getVal(fg, :x1))

setfreeze!(fg, [:x3, :x5])
tree, smtasks, hists = solveTree!(fg, tree)#, recordcliqs=ls(fg));
# csmAnimate(tree, hists, frames=1)

@test lm0 == getVal(fg, :lm0) #Still Frozen
@test X1 != getVal(fg, :x1) #not frozen

for i = [2,4,6]
    @test areCliqVariablesAllMarginalized(fg, tree.cliques[i])
end

end