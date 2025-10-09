using Test
using IncrementalInference

##

@testset "Testing basic marginalization" begin

## linear 6 
N=6
fg = generateGraph_LineStep(N; 
                            graphinit=false,
                            poseEvery=1, 
                            landmarkEvery=N+1, 
                            posePriorsAt=[0],
                            landmarkPriorsAt=[], 
                            sightDistance=N+1)

deleteFactor!.(fg, [Symbol("x$(i)lm0f1") for i=1:(N-1)])

# tree = buildTreeReset!(fg, show=true, drawpdf=true)
tree = solveTree!(fg)

for var in sortDFG(ls(fg))
    sppe = calcMeanMaxSuggested(fg, var, :default).suggested
    println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.35)
end


defaultFixedLagOnTree!(fg, 6)
# Back up data from these two poses so we can compare them once we solve again.
lm0 = deepcopy(getVal(fg, :lm0))
X0 = deepcopy(getVal(fg, :x0))
X1 = deepcopy(getVal(fg, :x1))

fifoFreeze!(fg)

##

tree = solveTree!(fg; recordcliqs=ls(fg));
for var in sortDFG(ls(fg))
    sppe = calcMeanMaxSuggested(fg, var, :default).suggested
    println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.35)
end

# although variables are marginalized, no cliques are
@test calcCliquesRecycled(tree) == (6,0,0,0)

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

smtasks = Task[];

eliminationOrder = [:x1, :x3, :x9, :x7, :x5, :lm0, :x8, :x4, :x2, :x6]

tree = solveTree!(fg; recordcliqs=ls(fg), smtasks, eliminationOrder);
hists = fetchCliqHistoryAll!(smtasks)

#clique 7 should be marginalized and therefor not do up or downsolve
@test calcCliquesRecycled(tree) == (7,1,0,0)
@test !(IIF.solveDown_StateMachine in getindex.(hists[7], 3))
@test !(IIF.solveUp_StateMachine in getindex.(hists[7], 3))

for var in sortDFG(ls(fg))
    sppe = calcMeanMaxSuggested(fg, var, :default).suggested
    println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.35)
end

@test lm0 == getVal(fg, :lm0) #Still Frozen
@test X1cmp == getVal(fg, :x1) #also now Frozen

unfreezeVariablesAll!(fg)

# freeze again
defaultFixedLagOnTree!(fg, 9)

tree = solveTree!(fg; eliminationOrder)

@test lm0 == getVal(fg, :lm0) #Still Frozen
@test X1cmp != getVal(fg, :x1) #not frozen

# freeze 6,8 to all marginalize clique 2
setfreeze!(fg, [:x6, :x8])
smtasks = Task[];
tree = solveTree!(fg; recordcliqs=ls(fg), smtasks, eliminationOrder);

hists = fetchCliqHistoryAll!(smtasks)

#clique 2 should be marginalized and therefor not do up or downsolve
@test calcCliquesRecycled(tree) == (7,1,0,0)
@test !(IIF.solveDown_StateMachine in getindex.(hists[2], 3))
@test !(IIF.solveUp_StateMachine in getindex.(hists[2], 3))
@test areCliqVariablesAllMarginalized(fg, getClique(tree,2))

tree = solveTree!(fg, tree; recordcliqs=ls(fg), eliminationOrder);
for var in sortDFG(ls(fg))
    sppe = calcMeanMaxSuggested(fg, var, :default).suggested
    println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.355)
end

@test calcCliquesRecycled(tree) == (7,1,6,0)

X1 = deepcopy(getVal(fg, :x1))

# to freeze clique 2,3,4
setfreeze!(fg, [:x4, :x5, :x7])

tree = solveTree!(fg, tree; recordcliqs=ls(fg), eliminationOrder);
# csmAnimate(tree, hists, frames=1)

@test calcCliquesRecycled(tree) == (7,3,4,0)

@test lm0 == getVal(fg, :lm0) #Still Frozen
@test X1 != getVal(fg, :x1) #not frozen

for i = [2,3,4]
    @test areCliqVariablesAllMarginalized(fg, getClique(tree, i))
end

end


@testset "Testing basic incremental recycle" begin

fg = generateGraph_LineStep(3; 
                                   poseEvery=1, 
                                   landmarkEvery=3, 
                                   posePriorsAt=[],
                                   landmarkPriorsAt=[0], 
                                   sightDistance=2,
                                   solverParams=SolverParams(algorithms=[:default, :parametric]))

getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true
# getSolverParams(fg).dbg = true

# tree = buildTreeReset!(fg, drawpdf=true, show=true)

eliminationOrder = [:lm3, :x0, :x3, :x1, :x2, :lm0]
tree = solveTree!(fg; recordcliqs=ls(fg), eliminationOrder); # , smtasks=smtasks

addFactor!(fg, [:lm3], Prior(Normal(3, 0.1)), graphinit=false)
smtasks = Task[]
tree = solveTree!(fg, tree; smtasks, recordcliqs=ls(fg), eliminationOrder);
hists = fetchCliqHistoryAll!(smtasks)

@test !(IIF.solveUp_StateMachine in getindex.(hists[3], 3))

for var in sortDFG(ls(fg))
    sppe = calcMeanMaxSuggested(fg, var, :default).suggested
    println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.35)
end


end


@testset "Testing incremental hex" begin
    
N=6
sfg = generateGraph_LineStep(N; 
                                   graphinit=false,
                                   poseEvery=1, 
                                   landmarkEvery=N+1, 
                                   posePriorsAt=[0],
                                   landmarkPriorsAt=[], 
                                   sightDistance=N+1,
                                   solverParams=SolverParams(algorithms=[:default, :parametric]))

deleteFactor!.(sfg, [Symbol("x$(i)lm0f1") for i=1:(N-1)])

vsyms = sortDFG(ls(sfg))
fsyms = sortDFG(lsf(sfg))

fg = deepcopyGraph(LocalDFG, sfg, vsyms[1:3], fsyms[1:3])

getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true
getSolverParams(fg).dbg = true
getSolverParams(fg).useMsgLikelihoods = true

tree = buildTreeReset!(fg)#, drawpdf=true, show=true)

smtasks = Task[]
tree = solveTree!(fg, tree; smtasks=smtasks, recordcliqs=ls(fg));

for var in sortDFG(ls(fg))
    sppe = calcMeanMaxSuggested(fg, var, :default).suggested
    # println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.35)
end

# add some more 
deepcopyGraph!(fg, sfg, vsyms[4:6], fsyms[4:6])

tree = solveTree!(fg, tree; smtasks=smtasks, recordcliqs=ls(fg));

for var in sortDFG(ls(fg))
    sppe = calcMeanMaxSuggested(fg, var, :default).suggested
    # println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.35)
end

# add some more 
deepcopyGraph!(fg, sfg, vsyms[7:8], fsyms[7:8])

tree = solveTree!(fg, tree; smtasks=smtasks, recordcliqs=ls(fg));

for var in sortDFG(ls(fg))
    sppe = calcMeanMaxSuggested(fg, var, :default).suggested
    # println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.35)
end

# add some more,close the loop 
deepcopyGraph!(fg, sfg, Symbol[], [fsyms[9]])

tree = solveTree!(fg, tree; smtasks=smtasks, recordcliqs=ls(fg));

for var in sortDFG(ls(fg))
    sppe = calcMeanMaxSuggested(fg, var, :default).suggested
    println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.35)
end

# force a reshuffle
addFactor!(fg, [:x4], Prior(Normal(4.1,0.1)))
tree = solveTree!(fg, tree; smtasks=smtasks, recordcliqs=ls(fg));

for var in sortDFG(ls(fg))
    sppe = calcMeanMaxSuggested(fg, var, :default).suggested
    # println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.35)
end

# and another reuse
addFactor!(fg, [:x4], Prior(Normal(3.9,0.1)))

tree = solveTree!(fg, tree; smtasks=smtasks, recordcliqs=ls(fg));

for var in sortDFG(ls(fg))
    sppe = calcMeanMaxSuggested(fg, var, :default).suggested
    # println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.355)
end

# printCSMHistoryLogical(hists)
# all except clique 1 should go for UPRECYCLED directly
#     | 1   x4             | 2   x2             | 3   x5             | 4   x3             | 5   lm0            | 6   x1             
# ----+--------------------+--------------------+--------------------+--------------------+--------------------+--------------------
# 1   | 1   setClique NULL | 3   setClique NULL | 4   setClique NULL | 5   setClique NULL | 6   setClique NULL | 7   setClique NULL 
# 2   | 2   buildCliq NULL | 8   buildCliq UPRE | 9   buildCliq UPRE | 10  buildCliq UPRE | 11  buildCliq UPRE | 12  buildCliq UPRE 
# 3   | 13  waitForUp NULL | 14  waitForUp UPRE | 15  waitForUp UPRE | 17  waitForUp UPRE | 19  waitForUp UPRE | 21  waitForUp UPRE 
# 4   | 34  preUpSolv NULL | 31  preUpSolv UPRE | 16  preUpSolv UPRE | 18  preUpSolv UPRE | 20  preUpSolv UPRE | 22  preUpSolv UPRE 
# 5   | 35  solveUp_S NULL | 32  postUpSol UPRE | 23  postUpSol UPRE | 24  postUpSol UPRE | 25  postUpSol UPRE | 26  postUpSol UPRE 
# 6   | 36  postUpSol UPSO | 33  waitForDo UPRE | 27  waitForDo UPRE | 28  waitForDo UPRE | 29  waitForDo UPRE | 30  waitForDo UPRE 
# 7   | 37  waitForDo UPSO | 39  solveDown UPRE | 40  solveDown UPRE | 43  solveDown UPRE | 41  solveDown UPRE | 44  solveDown UPRE 
# 8   | 38  solveDown DOWN | 47  updateFro DOWN | 46  updateFro DOWN | 52  updateFro DOWN | 50  updateFro DOWN | 54  updateFro DOWN 
# 9   | 42  updateFro DOWN | 49  exitState DOWN | 48  exitState DOWN | 53  exitState DOWN | 51  exitState DOWN | 55  exitState DOWN 
# 10  | 45  exitState DOWN | 


end