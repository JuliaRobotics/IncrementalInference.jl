using Test

using DistributedFactorGraphs
using IncrementalInference
using LieGroups
using LinearAlgebra

@testset "Test consolidation of factors #467" begin
  fg = generateGraph_LineStep(20, poseEvery=1, landmarkEvery=4, posePriorsAt=collect(0:7), sightDistance=2, solverParams=SolverParams(algorithms=[:default, :parametric]))
  IIF.prepare!(fg, IIF.NLLSSolver(), :parametric)
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

IIF.prepare!(fg, IIF.NLLSSolver(), :parametric)
vardict, result, flatvars, Σ = IIF.solveConditionalsParametric(fg, [:x1])
v1 = vardict[:x1]
@test isapprox(v1.val, [1.1], atol=1e-3)
# TODO what should the covariance be, should covariance on :x0 not influence it?
@test isapprox(v1.cov, [1.44;;], atol=1e-3)

initVariable!(fg, :x1, Normal(v1.val[1], sqrt(v1.cov[1])), :parametric)

addVariable!(fg, :x2, ContinuousScalar)
addFactor!(fg, [:x0,:x2], LinearRelative(Normal(2.0, 0.5)))
addFactor!(fg, [:x1,:x2], LinearRelative(Normal(1.1, 0.5)))

IIF.prepare!(fg, IIF.NLLSSolver(), :parametric)
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
IIF.prepare!(fg, IIF.NLLSSolver(), :parametric)
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
@test length(DFG.refMeans(v0.states[:parametric])[1]) === 1
@test isapprox(DFG.refMeans(v0.states[:parametric])[1][1], 0.0, atol = 1e-4)

@test IIF.autoinitParametric!(fg, :x1)

v0 = getVariable(fg,:x1)
@test length(DFG.refMeans(v0.states[:parametric])[1]) === 1
@test isapprox(DFG.refMeans(v0.states[:parametric])[1][1], 1.0, atol = 1e-4)

initAll!(fg)
IIF.initParametricFrom!(fg)

#
v0 = getVariable(fg,:x0)
@test length(DFG.refMeans(v0.states[:parametric])[1]) === 1
@test isapprox(DFG.refMeans(v0.states[:parametric])[1][1], 0.0, atol = 0.1)
v1 = getVariable(fg,:x1)
@test isapprox(DFG.refMeans(v1.states[:parametric])[1][1], 1.0, atol = 0.1)

##

fg = generateGraph_LineStep(10, vardims=2, poseEvery=1, landmarkEvery=3, posePriorsAt=Int[0,5,10], sightDistance=3, solverParams=SolverParams(algorithms=[:default, :parametric]))
    # addFactor!(fg, [:x5; :x15], LinearRelative(Normal(10, 0.1)))
    # addFactor!(fg, [:x15; :x25], LinearRelative(Normal(10, 0.1)))

#to manually check all factors
# foreach(fct->println(fct.label, ": ", getObservation(fct).Z), getFactors(fg))

# @profiler d,st = IIF.solveGraphParametric(fg)
IIF.prepare!(fg, IIF.NLLSSolver(), :parametric)
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

foreach(x->DFG.refMeans(DFG.getState(fg, x.first, :parametric))[1] = x.second, pairs(d))


# getSolverParams(fg).dbg=true
# getSolverParams(fg).drawtree=true
# getSolverParams(fg).async = true
getSolverParams(fg).graphinit = false

tree2 = IIF.solveTree!(fg; algorithm = :parametric) #, recordcliqs=ls(fg))


for i in 0:10
  sym = Symbol("x",i)
  @show val = DFG.refMeans(DFG.getState(fg, sym, :parametric))
  @test isapprox(val[1][1], i, atol=1e-3)
  @test isapprox(val[1][2], i, atol=1e-3)
end

##

# Print answers
if false
vsds = DFG.getState.(getVariables(fg), :parametric)
foreach(v->println(v.label, ": ", DFG.refMeans(DFG.getState(v, :parametric))), sort!(getVariables(fg), by=getLabel, lt=natural_lt))
end


## #################################################################

fg = LocalDFG( solverParams=SolverParams(algorithms=[:default, :parametric]))
# fg = LocalDFG{SolverParams}( solverParams=SolverParams())
N = 100
getSolverParams(fg).N = N
graphinit = false

addVariable!(fg, :x0, ContinuousScalar, N=N) # autoinit = graphinit
addFactor!(fg, [:x0], Prior(Normal(-1.0, 1.0)))

addVariable!(fg, :x1, ContinuousScalar, N=N) # autoinit = graphinit

addVariable!(fg, :x2, ContinuousScalar, N=N) # autoinit = graphinit
addFactor!(fg, [:x2], Prior(Normal(+1.0, 1.0)))

addFactor!(fg, [:x0; :x1], LinearRelative(Normal(0.0, 1e-1)), graphinit=graphinit)
addFactor!(fg, [:x1; :x2], LinearRelative(Normal(0.0, 1e-1)), graphinit=graphinit)




foreach(fct->println(fct.label, ": ", getObservation(fct).Z), getFactors(fg))
IIF.prepare!(fg, IIF.NLLSSolver(), :parametric)
M, labels, minimizer, Σ = IIF.solveGraphParametric(fg)
d = Dict(labels.=>minimizer)

foreach(println, d)
@test isapprox(d[:x0][1][1], -0.01, atol=1e-3)
@test isapprox(d[:x1][1][1], 0.0, atol=1e-3)
@test isapprox(d[:x2][1][1], 0.01, atol=1e-3)


##

foreach(x->DFG.refMeans(DFG.getState(getVariable(fg,x.first),:parametric))[1] = x.second, pairs(d))

# task = @async begin
  #   global tree2
  #   global smt
  #   global hist
#force message passing with manual variable order
tree2 = solveTree!(fg; algorithm=:parametric, eliminationOrder=[:x0, :x2, :x1])
# end
foreach(v->println(v.label, ": ", DFG.refMeans(DFG.getState(v, :parametric))), getVariables(fg))

@test isapprox(DFG.refMeans(getVariable(fg,:x0).states[:parametric])[1][1], -0.01, atol=1e-3)
@test isapprox(DFG.refMeans(getVariable(fg,:x1).states[:parametric])[1][1], 0.0, atol=1e-3)
@test isapprox(DFG.refMeans(getVariable(fg,:x2).states[:parametric])[1][1], 0.01, atol=1e-3)

## ##############################################################################
## multiple sections

fg = generateGraph_LineStep(10, poseEvery=1, landmarkEvery=10, posePriorsAt=Int[0,10], sightDistance=5, solverParams=SolverParams(algorithms=[:default, :parametric]))
# break fg in 2
deleteFactor!(fg, :x5x6f1)
# plotDFG(fg)

#check all factors
# foreach(fct->println(fct.label, ": ", getObservation(fct).Z), getFactors(fg))

# @profiler d,st = IIF.solveGraphParametric(fg)
IIF.prepare!(fg, IIF.NLLSSolver(), :parametric)
M, labels, minimizer, Σ = IIF.solveGraphParametric(fg)
d = Dict(labels.=>minimizer)
if false
foreach(println, d)
end
for i in 0:10
  sym = Symbol("x",i)
  @test isapprox(d[sym][1], i, atol=1e-6)
end

foreach(x->DFG.refMeans(DFG.getState(getVariable(fg,x.first),:parametric))[1] = x.second, pairs(d))

getSolverParams(fg).graphinit = false
tree2 = IIF.solveTree!(fg; algorithm=:parametric)

# print results
if false
vsds = DFG.getState.(getVariables(fg), :parametric)
foreach(v->println(v.label, ": ", DFG.refMeans(DFG.getState(v, :parametric))), getVariables(fg))
end

for i in 0:10
  sym = Symbol("x",i)
  var = getVariable(fg,sym)
  val = DFG.refMeans(var.states[:parametric])
  #TODO investigate why tolarance degraded (its tree related and not bad enough to worry now)
  @test isapprox(val[1][1], i, atol=5e-4) 
end

##

end


@testset "initAll!(fg, :parametric)" begin
##

fg = generateGraph_LineStep(7, poseEvery=1, landmarkEvery=0, posePriorsAt=collect(0:7), sightDistance=2, solverParams=SolverParams(graphinit=false), graphinit=false)

IIF.prepare!(fg, IIF.NLLSSolver(), :parametric)
@test (l->!isInitialized(fg, l, :parametric)).(ls(fg)) |> all

initAll!(fg, :parametric)

@test (l->isInitialized(fg, l, :parametric)).(ls(fg)) |> all


##
end

# test/testParametricUninitializable.jl
## Define a ternary "biased relative" factor: x_j = x_i + z + b
# The bias `b` is only observable through these factors (no direct prior).
struct BiasedLinearRelative{T <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
    Z::T
end

DFG.getManifold(::IIF.InstanceType{BiasedLinearRelative}) = LieGroups.TranslationGroup(1)

# residual: z - (x2 - x1 - b)
function (cf::CalcFactor{<:BiasedLinearRelative})(z, x1, x2, b)
    return z .- (x2 .- x1 .- b)
end

##
@testset "Parametric: uninitializable variable (ternary bias factor)" begin
    fg = initfg()
    fg.solverParams.graphinit = false

    # Chain: x0 --[biased]--> x1 --[biased]--> x2
    # with shared bias variable :b (no prior on :b)
    addVariable!(fg, :x0, ContinuousScalar)
    addVariable!(fg, :x1, ContinuousScalar)
    addVariable!(fg, :x2, ContinuousScalar)
    addVariable!(fg, :b, ContinuousScalar)

    addFactor!(fg, [:x0], Prior(Normal(0.0, 0.1)))
    addFactor!(fg, [:x2], Prior(Normal(2.5, 0.1)))

    # Biased relative factors: x1 = x0 + 1.0 + b, x2 = x1 + 1.0 + b
    # True solution: 1.0+b → 2*(1.0+b)=2.5 → b=0.25
    addFactor!(fg, [:x0, :x1, :b], BiasedLinearRelative(Normal(1.0, 0.1)))
    addFactor!(fg, [:x1, :x2, :b], BiasedLinearRelative(Normal(1.0, 0.1)))

    # autoinitParametric! should NOT crash on :b even though it's locally under-constrained
    IIF.autoinitParametric!(fg)

    # The global parametric solve should still work and find the correct solution
    M, v, r, Λ = IIF.solveGraphParametric!(fg; init=false)

    x0 = DFG.refMeans(getState(fg, :x0, :parametric))[1]
    x1 = DFG.refMeans(getState(fg, :x1, :parametric))[1]
    x2 = DFG.refMeans(getState(fg, :x2, :parametric))[1]
    b  = DFG.refMeans(getState(fg, :b, :parametric))[1]

    @test isapprox(x0[1], 0.0, atol=0.05)
    @test isapprox(x2[1], 2.5, atol=0.05)
    @test isapprox(b[1], 0.25, atol=0.05)
    @test isapprox(x1[1], x0[1] + 1.0 + b[1], atol=0.05)
end

"""
    PartialExpCoordPrior

A partial prior that constrains specific exponential coordinates of a Lie group variable.

Mathematically, this factor applies a prior to a subset of the tangent coordinates `vee(log(G, g))`. 
It acts as a locally valid submersion on any Lie group, provided the variable remains within 
the injectivity radius where the parameterization in exponential coordinates is well-defined.
The `partial` tuple selects which coordinates of the vee representation are observed.
"""
struct PartialExpCoordPrior{G <: LieGroups.AbstractLieGroup, T <: IIF.SamplableBelief, P <: Tuple} <: IIF.AbstractPriorObservation
    G::G
    Z::T
    partial::P
end

# Factor manifold is the residual space: ℝ^k where k = length(partial)
DFG.getManifold(pp::PartialExpCoordPrior) = LieGroups.TranslationGroup(length(pp.partial))

function (cf::CalcFactor{<:PartialExpCoordPrior})(z, x1)
    G = cf.factor.G
    # Get exponential coordinates
    X = log(G, x1) 
    Xc = vee(LieAlgebra(G), X)
    return z .- Xc[collect(cf.factor.partial)]   # Residual on selected coords
end

@testset "Parametric: PartialExpCoordPrior on 2D variable (locally rank-deficient)" begin
    fg = initfg()
    fg.solverParams.graphinit = false

    G = LieGroups.TranslationGroup(2)

    # x0 has partial prior on x-coord only (y unconstrained locally)
    # x2 has partial prior on y-coord only (x unconstrained locally)
    # LinearRelative{2} chain makes the full graph solvable
    addVariable!(fg, :x0, ContinuousEuclid{2})
    addVariable!(fg, :x1, ContinuousEuclid{2})
    addVariable!(fg, :x2, ContinuousEuclid{2})

    # x0: only x-coordinate known via partial prior on coord 1
    addFactor!(fg, [:x0], PartialExpCoordPrior(G, Normal(0.0, 0.1), (1,)))
    # x2: only y-coordinate known via partial prior on coord 2
    addFactor!(fg, [:x2], PartialExpCoordPrior(G, Normal(3.0, 0.1), (2,)))

    # Relative factors that constrain both dimensions
    addFactor!(fg, [:x0, :x1], LinearRelative{2}(MvNormal([1.0, 1.0], 0.1*I(2))))
    addFactor!(fg, [:x1, :x2], LinearRelative{2}(MvNormal([1.0, 1.0], 0.1*I(2))))

    IIF.autoinitParametric!(fg)

    M, v, r, Λ = IIF.solveGraphParametric!(fg; init=false)

    x0 = DFG.refMeans(getState(fg, :x0, :parametric))[1]
    x1 = DFG.refMeans(getState(fg, :x1, :parametric))[1]
    x2 = DFG.refMeans(getState(fg, :x2, :parametric))[1]

    # x0[1] ≈ 0.0 (from prior), x2[2] ≈ 3.0 (from prior)
    # Propagation: x0[2] = x2[2] - 2.0 = 1.0, x2[1] = x0[1] + 2.0 = 2.0
    @test isapprox(x0, [0.0, 1.0], atol=0.05)
    @test isapprox(x1, [1.0, 2.0], atol=0.05)
    @test isapprox(x2, [2.0, 3.0], atol=0.05)
end
