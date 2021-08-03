using DistributedFactorGraphs
using IncrementalInference
using Manifolds
using StaticArrays
using Test

##

@testset "Test SpecialOrthogonal(2) prior" begin

##

Base.convert(::Type{<:Tuple}, M::SpecialOrthogonal{2}) = (:Circular,)
Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{SpecialOrthogonal{2}})  = (:Circular,)

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

##

doautoinit!(fg, :x0)

vnd = getVariableSolverData(fg, :x0)
@test all(isapprox.(mean(vnd.val), [1 0; 0 1], atol=0.1))
@test all(is_point.(Ref(M), vnd.val))


##
v1 = addVariable!(fg, :x1, SpecialOrthogonal2)
mf = ManifoldFactor(SpecialOrthogonal(2), MvNormal([pi], [0.01]))
f = addFactor!(fg, [:x0, :x1], mf)

doautoinit!(fg, :x1)

##
# smtasks = Task[]
solveTree!(fg) #; smtasks, verbose=true, recordcliqs=ls(fg))
# hists = fetchCliqHistoryAll!(smtasks);

##

end


@testset "Test SpecialOrthogonal(3) prior" begin

##

Base.convert(::Type{<:Tuple}, M::SpecialOrthogonal{3}) = (:Euclid, :Euclid, :Euclid)
Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{SpecialOrthogonal{3}})  =  (:Euclid, :Euclid, :Euclid)

# @defVariable SO3 SpecialOrthogonal(3) @MMatrix([1.0 0.0; 0.0 1.0])
@defVariable SO3 SpecialOrthogonal(3) [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

M = getManifold(SO3)
@test M == SpecialOrthogonal(3)
pT = getPointType(SO3)
# @test pT == MMatrix{2, 2, Float64, 4}
@test pT == Matrix{Float64}
pϵ = getPointIdentity(SO3)
@test pϵ == [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

@test is_point(getManifold(SO3), getPointIdentity(SO3))

fg = initfg()

v0 = addVariable!(fg, :x0, SO3)

mp = ManifoldPrior(SpecialOrthogonal(3), SA[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], MvNormal([0.01, 0.01, 0.01]))
p = addFactor!(fg, [:x0], mp)

doautoinit!(fg, :x0)

vnd = getVariableSolverData(fg, :x0)
@test all(isapprox.( mean(SpecialOrthogonal(3),vnd.val), [1 0 0; 0 1 0; 0 0 1], atol=0.01))
@test all(is_point.(Ref(M), vnd.val))

points = sampleFactor(fg, :x0f1, 100)[1]
IIF.calcCovarianceBasic(SpecialOrthogonal(3), points)

##
v1 = addVariable!(fg, :x1, SO3)
mf = ManifoldFactor(SpecialOrthogonal(3), MvNormal([0.01,0.01,0.01], [0.01,0.01,0.01]))
f = addFactor!(fg, [:x0, :x1], mf)

doautoinit!(fg, :x1)

vnd = getVariableSolverData(fg, :x1)
@test all(isapprox.( mean(SpecialOrthogonal(3),vnd.val), [0.9999 -0.00995 0.01005; 0.01005 0.9999 -0.00995; -0.00995 0.01005 0.9999], atol=0.01))
@test all(is_point.(Ref(M), vnd.val))

# smtasks = Task[]
solveTree!(fg) # ; smtasks, verbose=true, recordcliqs=ls(fg))

# test them again after solve
vnd = getVariableSolverData(fg, :x0)
@test all(isapprox.( mean(SpecialOrthogonal(3),vnd.val), [1 0 0; 0 1 0; 0 0 1], atol=0.01))
@test all(is_point.(Ref(M), vnd.val))

vnd = getVariableSolverData(fg, :x1)
@test all(isapprox.( mean(SpecialOrthogonal(3),vnd.val), [0.9999 -0.00995 0.01005; 0.01005 0.9999 -0.00995; -0.00995 0.01005 0.9999], atol=0.01))
@test all(is_point.(Ref(M), vnd.val))

##
end