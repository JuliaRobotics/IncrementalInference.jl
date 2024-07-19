using Test

using DistributedFactorGraphs
using IncrementalInference

##

@testset "Test consolidation of factors #467" begin
  fg = generateGraph_LineStep(20, poseEvery=1, landmarkEvery=4, posePriorsAt=collect(0:7), sightDistance=2, solverParams=SolverParams(algorithms=[:default, :parametric]))

  M, labels, minimizer, Σ = IIF.solveGraphParametric(fg)
  d = Dict(labels.=>minimizer)
  for i in 0:20
    sym = Symbol("x",i)
    @test isapprox(d[sym][1], i, atol=1e-6)
  end
  
  for i in 0:4:20
    sym = Symbol("lm",i)
    @test isapprox(d[sym][1], i, atol=1e-6)
  end
  
end

##
@testset "Parametric Tests" begin
fg = LocalDFG(solverParams=SolverParams(algorithms=[:default, :parametric]))

addVariable!(fg, :x0, ContinuousScalar)
initVariable!(fg, :x0, Normal(0.1,1.1), :parametric)

addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0,:x1], LinearRelative(Normal(1.0, 1.2)))

vardict, result, flatvars, Σ = IIF.solveConditionalsParametric(fg, [:x1])
v1 = vardict[:x1]
@test isapprox(v1.val, [1.1], atol=1e-3)
# TODO what should the covariance be, should covariance on :x0 not influence it?
@test isapprox(v1.cov, [1.44;;], atol=1e-3)

initVariable!(fg, :x1, Normal(v1.val[1], sqrt(v1.cov[1])), :parametric)

addVariable!(fg, :x2, ContinuousScalar)
addFactor!(fg, [:x0,:x2], LinearRelative(Normal(2.0, 0.5)))
addFactor!(fg, [:x1,:x2], LinearRelative(Normal(1.1, 0.5)))


vardict, result, flatvars, Σ = IIF.solveConditionalsParametric(fg, [:x2])
v2 = vardict[:x2]

@test isapprox(v2.val, [2.15], atol=1e-3)
# TODO what should the covariance be?
@test isapprox(v2.cov, [0.125;;], atol=1e-3)
initVariable!(fg, :x2, Normal(v2.val[1], sqrt(v2.cov[1])), :parametric)

addFactor!(fg, [:x0], Prior(Normal(0.1,1.1)))
IIF.solveGraphParametric!(fg; is_sparse=false)

end

@testset "Parametric Tests" begin

##
fg = generateGraph_LineStep(7, poseEvery=1, landmarkEvery=0, posePriorsAt=collect(0:7), sightDistance=2, solverParams=SolverParams(algorithms=[:default, :parametric]))

M, labels, minimizer, Σ = IIF.solveGraphParametric(fg)
d = Dict(labels.=>minimizer)

for i in 0:7
  sym = Symbol("x",i)
  @test isapprox(d[sym][1], i, atol=1e-6)
end


##

fg = generateGraph_LineStep(2, graphinit=true, vardims=1, poseEvery=1, landmarkEvery=0, posePriorsAt=Int[0], sightDistance=3, solverParams=SolverParams(algorithms=[:default, :parametric]))

@test IIF.autoinitParametric!(fg, :x0)

v0 = getVariable(fg,:x0)
@test length(v0.solverDataDict[:parametric].val[1]) === 1
@test isapprox(v0.solverDataDict[:parametric].val[1][1], 0.0, atol = 1e-4)

@test IIF.autoinitParametric!(fg, :x1)

v0 = getVariable(fg,:x1)
@test length(v0.solverDataDict[:parametric].val[1]) === 1
@test isapprox(v0.solverDataDict[:parametric].val[1][1], 1.0, atol = 1e-4)


IIF.initParametricFrom!(fg)

#
v0 = getVariable(fg,:x0)
@test length(v0.solverDataDict[:parametric].val[1]) === 1
@test isapprox(v0.solverDataDict[:parametric].val[1][1], 0.0, atol = 0.1)
v1 = getVariable(fg,:x1)
@test isapprox(v1.solverDataDict[:parametric].val[1][1], 1.0, atol = 0.1)

##

fg = generateGraph_LineStep(10, vardims=2, poseEvery=1, landmarkEvery=3, posePriorsAt=Int[0,5,10], sightDistance=3, solverParams=SolverParams(algorithms=[:default, :parametric]))
    # addFactor!(fg, [:x5; :x15], LinearRelative(Normal(10, 0.1)))
    # addFactor!(fg, [:x15; :x25], LinearRelative(Normal(10, 0.1)))

#to manually check all factors
# foreach(fct->println(fct.label, ": ", getFactorType(fct).Z), getFactors(fg))

# @profiler d,st = IIF.solveGraphParametric(fg)
M, labels, minimizer, Σ = IIF.solveGraphParametric(fg)
d = Dict(labels.=>minimizer)

for i in 0:10
  sym = Symbol("x",i)
  @test isapprox(d[sym][1], i, atol=1e-6)
  @test isapprox(d[sym][2], i, atol=1e-6)
end

# print results out
if false
  foreach(println, d)
end


##

foreach(x->getSolverData(getVariable(fg,x.first),:parametric).val[1] = x.second, pairs(d))


# getSolverParams(fg).dbg=true
# getSolverParams(fg).drawtree=true
# getSolverParams(fg).async = true
getSolverParams(fg).graphinit = false

tree2 = IIF.solveTree!(fg; algorithm = :parametric) #, recordcliqs=ls(fg))


for i in 0:10
  sym = Symbol("x",i)
  var = getVariable(fg,sym)
  @show val = var.solverDataDict[:parametric].val
  @test isapprox(val[1][1], i, atol=1e-3)
  @test isapprox(val[1][2], i, atol=1e-3)
end

##

# Print answers
if false
vsds = DFG.getSolverData.(getVariables(fg), :parametric)
foreach(v->println(v.label, ": ", DFG.getSolverData(v, :parametric).val), sort!(getVariables(fg), by=getLabel, lt=natural_lt))
end


## #################################################################

fg = LocalDFG( solverParams=SolverParams(algorithms=[:default, :parametric]))
# fg = LocalDFG{SolverParams}( solverParams=SolverParams())
N = 100
fg.solverParams.N = N
graphinit = false

addVariable!(fg, :x0, ContinuousScalar, N=N) # autoinit = graphinit
addFactor!(fg, [:x0], Prior(Normal(-1.0, 1.0)))

addVariable!(fg, :x1, ContinuousScalar, N=N) # autoinit = graphinit

addVariable!(fg, :x2, ContinuousScalar, N=N) # autoinit = graphinit
addFactor!(fg, [:x2], Prior(Normal(+1.0, 1.0)))

addFactor!(fg, [:x0; :x1], LinearRelative(Normal(0.0, 1e-1)), graphinit=graphinit)
addFactor!(fg, [:x1; :x2], LinearRelative(Normal(0.0, 1e-1)), graphinit=graphinit)




foreach(fct->println(fct.label, ": ", getFactorType(fct).Z), getFactors(fg))

M, labels, minimizer, Σ = IIF.solveGraphParametric(fg)
d = Dict(labels.=>minimizer)

foreach(println, d)
@test isapprox(d[:x0][1][1], -0.01, atol=1e-3)
@test isapprox(d[:x1][1][1], 0.0, atol=1e-3)
@test isapprox(d[:x2][1][1], 0.01, atol=1e-3)


##

foreach(x->getSolverData(getVariable(fg,x.first),:parametric).val[1] = x.second, pairs(d))

# fg.solverParams.showtree = true
# fg.solverParams.drawtree = true
# fg.solverParams.dbg = true
# fg.solverParams.graphinit = false

# task = @async begin
  #   global tree2
  #   global smt
  #   global hist
#force message passing with manual variable order
tree2 = solveTree!(fg; algorithm=:parametric, eliminationOrder=[:x0, :x2, :x1])
# end
foreach(v->println(v.label, ": ", DFG.getSolverData(v, :parametric).val), getVariables(fg))

@test isapprox(getVariable(fg,:x0).solverDataDict[:parametric].val[1][1], -0.01, atol=1e-3)
@test isapprox(getVariable(fg,:x1).solverDataDict[:parametric].val[1][1], 0.0, atol=1e-3)
@test isapprox(getVariable(fg,:x2).solverDataDict[:parametric].val[1][1], 0.01, atol=1e-3)

## ##############################################################################
## multiple sections

fg = generateGraph_LineStep(10, poseEvery=1, landmarkEvery=10, posePriorsAt=Int[0,10], sightDistance=5, solverParams=SolverParams(algorithms=[:default, :parametric]))
# break fg in 2
deleteFactor!(fg, :x5x6f1)
# plotDFG(fg)

#check all factors
# foreach(fct->println(fct.label, ": ", getFactorType(fct).Z), getFactors(fg))

# @profiler d,st = IIF.solveGraphParametric(fg)
M, labels, minimizer, Σ = IIF.solveGraphParametric(fg)
d = Dict(labels.=>minimizer)
if false
foreach(println, d)
end
for i in 0:10
  sym = Symbol("x",i)
  @test isapprox(d[sym][1], i, atol=1e-6)
end

foreach(x->getSolverData(getVariable(fg,x.first),:parametric).val[1] = x.second, pairs(d))

# fg.solverParams.showtree = true
# fg.solverParams.drawtree = true
# fg.solverParams.dbg = false
getSolverParams(fg).graphinit = false
tree2 = IIF.solveTree!(fg; algorithm=:parametric)

# print results
if false
vsds = DFG.getSolverData.(getVariables(fg), :parametric)
foreach(v->println(v.label, ": ", DFG.getSolverData(v, :parametric).val), getVariables(fg))
end

for i in 0:10
  sym = Symbol("x",i)
  var = getVariable(fg,sym)
  val = var.solverDataDict[:parametric].val
  #TODO investigate why tolarance degraded (its tree related and not bad enough to worry now)
  @test isapprox(val[1][1], i, atol=5e-4) 
end

##

end


@testset "initAll!(fg, :parametric)" begin
##

fg = generateGraph_LineStep(7, poseEvery=1, landmarkEvery=0, posePriorsAt=collect(0:7), sightDistance=2, solverParams=SolverParams(graphinit=false), graphinit=false)

@test (l->!isInitialized(fg, l, :parametric)).(ls(fg)) |> all

initAll!(fg, :parametric)

@test (l->isInitialized(fg, l, :parametric)).(ls(fg)) |> all


##
end


#
