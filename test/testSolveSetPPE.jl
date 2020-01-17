# test that PPE values are update with solve, see issue #548

using Test
using IncrementalInference
using DistributedFactorGraphs

@testset "test PPE update during solve" begin

fg = generateCanonicalFG_Kaess(graphinit=true)
ensureAllInitialized!(fg)

# list of variables to check
vars = getVariableIds(fg)

# fetch values before solve
before = Dict()
for vs in vars
  before[vs] = getVariablePPE(getVariable(fg, vs)) |> getSuggestedPPE
end

# do the solve
# getSolverParams(fg).dbg = true

# tree = wipeBuildNewTree!(fg)
# drawTree(tree, show=true)

# solveCliq!(fg, tree, :l2)
# solveCliq!(fg, tree, :x3)
# solveCliq!(fg, tree, :x2)

solveTree!(fg)


after = Dict()
for vs in vars
  after[vs] = getVariablePPE(getVariable(fg, vs)) |> getSuggestedPPE
end

# before and after should be noticably different, because first inferred values have been found
for vs in vars
  errd = norm(before[vs] - after[vs])
  # @show vs, errd
  @test 1e-5 < errd
end

# force recalc and update each PPE
force = Dict()
for vs in vars
  setVariablePosteriorEstimates!(fg, vs)
  force[vs] = getVariablePPE(getVariable(fg, vs)) |> getSuggestedPPE
  # these need to be close to the same as after
  errd = norm(force[vs] - after[vs])
  # @show vs, errd
  @test errd < 0.1
end



## suspect cliqSubFg updated, but not back to main dfg object... test via load graph


end
