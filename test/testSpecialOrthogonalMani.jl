using DistributedFactorGraphs
using IncrementalInference
using Manifolds
using StaticArrays
using Test

# @testset "Test SpecialOrthogonal(2) prior" begin

Base.convert(::Type{<:Tuple}, M::SpecialOrthogonal{2}) = (:Euclid,)
Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{SpecialOrthogonal{2}})  = (:Euclid,)

# @defVariable SpecialOrthogonal2 SpecialOrthogonal(2) @MMatrix([1.0 0.0; 0.0 1.0])
@defVariable SpecialOrthogonal2 SpecialOrthogonal(2) [1.0 0.0; 0.0 1.0]

M = getManifold(SpecialOrthogonal2)
@test M == SpecialOrthogonal(2)
pT = getPointType(SpecialOrthogonal2)
# @test pT == MMatrix{2, 2, Float64, 4}
@test pT == Matrix{Float64}
pϵ = getPointIdentity(SpecialOrthogonal2)
@test pϵ == [1.0 0.0; 0.0 1.0]

@test is_point(getManifold(SpecialOrthogonal2), getPointIdentity(SpecialOrthogonal2))

fg = initfg()

v0 = addVariable!(fg, :x0, SpecialOrthogonal2)

mp = ManifoldPrior(SpecialOrthogonal(2), SA[1.0 0.0; 0.0 1.0], MvNormal([0.01]))
p = addFactor!(fg, [:x0], mp)

doautoinit!(fg, :x0)

vnd = getVariableSolverData(fg, :x0)
@test all(isapprox.(mean(vnd.val), [1 0; 0 1], atol=0.1))
@test all(is_point.(Ref(M), vnd.val))


##
v1 = addVariable!(fg, :x1, SpecialOrthogonal2)
mf = ManifoldFactor(SpecialOrthogonal(2), MvNormal([pi], [0.01]))
f = addFactor!(fg, [:x0, :x1], mf)

@test_broken doautoinit!(fg, :x1)

##
# Debugging SpecialOrthogonal error
smtasks = Task[]
@test_broken solveTree!(fg; smtasks, verbose=true, recordcliqs=ls(fg))
# hists = fetchCliqHistoryAll!(smtasks);

end