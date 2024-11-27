# using Revise

using DistributedFactorGraphs
using IncrementalInference
using Test

##

# @testset "deprecation of old api" begin

# LinearRelative(Normal(2.0, 0.1))

# end


@testset "testing compare functions for variables and factors..." begin

##

fg = initfg()

v1 = addVariable!(fg, :x0, ContinuousScalar)
f1 = addFactor!(fg, [:x0;], Prior(Normal()))

@test compareVariable(v1,v1)
@test compareFactor(f1,f1)

v2 = addVariable!(fg, :x1, ContinuousScalar)
f2 = addFactor!(fg, [:x0;:x1], LinearRelative(Normal(2.0, 0.1)))

fg2 = deepcopy(fg)

@test !compareVariable(v1,v2)

# not testing different factors in this way
# @test !compareFactor(f1,f2)

@test compareAllVariables(fg, fg)

@test compareAllVariables(fg, fg2)

@test compareSimilarVariables(fg, fg)
@test compareSimilarVariables(fg, fg2)

@test compareSimilarFactors(fg, fg)
@test compareSimilarFactors(fg, fg2; skip=[:particleidx])

@test compareFactorGraphs(fg, fg)
@test compareFactorGraphs(fg, fg2; skip=[:particleidx; :varidx])

# easier error messages
getSolverParams(fg).multiproc = false

tree = solveTree!(fg)

x1a = getVariable(fg, :x0)
x1b = getVariable(fg2, :x0)

@test !compareVariable(x1a, x1b, skipsamples=false)

@test !compareSimilarVariables(fg, fg2, skipsamples=false)
@test !compareSimilarFactors(fg, fg2, skipsamples=false, skip=[:measurement;])

@test compareFactorGraphs(fg, fg)
@test !compareFactorGraphs(fg, fg2, skipsamples=false)

initAll!(fg2)

@test compareSimilarVariables(fg, fg2, skipsamples=true, skip=Symbol[:initialized;:infoPerCoord;:ppeDict;:solvedCount])
# fg2 has been solved, so it should fail on the estimate dictionary
@test !compareSimilarVariables(fg, fg2, skipsamples=true, skip=Symbol[:initialized;:infoPerCoord])

tree = buildTreeReset!(fg2)

# Expect ccw to reflect different numerics since fg and fg2 have different numeric solutions
Al = IIF._getCCW(fg,  getLabel(f2))
Bl = IIF._getCCW(fg2, getLabel(f2))
field = :varValsAll
@test !compareField(Al, Bl, field)

@test compareSimilarFactors(fg, fg2, skipsamples=true, skipcompute=true, skip=[:fullvariables; :varValsAll; :particleidx])

@test !compareSimilarFactors(fg, fg2, skipsamples=true, skipcompute=false)


##

end

@testset "test subgraph functions..." begin

##

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal()))

addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0;:x1], LinearRelative(Normal(2.0, 0.1)))

addVariable!(fg, :x2, ContinuousScalar)
addFactor!(fg, [:x1;:x2], LinearRelative(Normal(4.0, 0.1)))

addVariable!(fg, :l1, ContinuousScalar)
addFactor!(fg, [:x1;:l1], LinearRelative(Rayleigh()))

sfg = buildSubgraph(fg, [:x0;:x1], 1) # distance=1 to include factors 

#FIXME JT - this doesn't make sense to pass, it is a subgraph so should it not rather be ⊂ [subset]?
# compareDFG(fg1, fg2, by=⊂, skip=...)
@test_broken compareFactorGraphs(fg, sfg, skip=[:labelDict;:addHistory;:logpath;:sessionLabel; :particleidx; :varidx])

# drawGraph(sfg)

##

end







#
