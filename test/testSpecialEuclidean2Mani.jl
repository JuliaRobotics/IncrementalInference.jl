using DistributedFactorGraphs
using IncrementalInference
using Manifolds
using StaticArrays
using Test

@testset "Test SpecialEuclidean(2)" begin

Base.convert(::Type{<:Tuple}, M::SpecialEuclidean{2}) = (:Euclid, :Euclid, :Circular)
Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{SpecialEuclidean{2}})  = (:Euclid, :Euclid, :Circular)

# @defVariable SpecialEuclidean2 SpecialEuclidean(2) ProductRepr(@MVector([0.0,0.0]), @MMatrix([1.0 0.0; 0.0 1.0]))
@defVariable SpecialEuclidean2 SpecialEuclidean(2) ProductRepr([0.0,0.0], [1.0 0.0; 0.0 1.0])

M = getManifold(SpecialEuclidean2)
@test M == SpecialEuclidean(2)
pT = getPointType(SpecialEuclidean2)
@test pT == ProductRepr{Tuple{Vector{Float64}, Matrix{Float64}}}
# @test pT == ProductRepr{Tuple{MVector{2, Float64}, MMatrix{2, 2, Float64, 4}}}
pϵ = getPointIdentity(SpecialEuclidean2)
# @test_broken pϵ == ProductRepr(@MVector([0.0,0.0]), @MMatrix([1.0 0.0; 0.0 1.0]))
@test all(isapprox.(pϵ,ProductRepr([0.0,0.0], [1.0 0.0; 0.0 1.0])).parts)

@test is_point(getManifold(SpecialEuclidean2), getPointIdentity(SpecialEuclidean2))

##
fg = initfg()

v0 = addVariable!(fg, :x0, SpecialEuclidean2)

# mp = ManifoldPrior(SpecialEuclidean(2), ProductRepr(@MVector([0.0,0.0]), @MMatrix([1.0 0.0; 0.0 1.0])), MvNormal([0.01, 0.01, 0.01]))
mp = ManifoldPrior(SpecialEuclidean(2), ProductRepr(@MVector([0.0,0.0]), @MMatrix([1.0 0.0; 0.0 1.0])), MvNormal([0.01, 0.01, 0.01]))
p = addFactor!(fg, [:x0], mp)

doautoinit!(fg, :x0)

##
vnd = getVariableSolverData(fg, :x0)
@test all(isapprox.(mean(vnd.val), ProductRepr([0.0,0.0], [1.0 0.0; 0.0 1.0]), atol=0.1).parts)
@test all(is_point.(Ref(M), vnd.val))

##
v1 = addVariable!(fg, :x1, SpecialEuclidean2)
mf = ManifoldFactor(SpecialEuclidean(2), MvNormal([1,2,pi/4], [0.01,0.01,0.01]))
f = addFactor!(fg, [:x0, :x1], mf)

doautoinit!(fg, :x1)

vnd = getVariableSolverData(fg, :x1)
@test all(isapprox.(mean(vnd.val), ProductRepr([1.0,2.0], [0.7071 -0.7071; 0.7071 0.7071]), atol=0.1).parts)
@test all(is_point.(Ref(M), vnd.val))

##
smtasks = Task[]
solveTree!(fg; smtasks, verbose=true, recordcliqs=ls(fg))
# hists = fetchCliqHistoryAll!(smtasks);

vnd = getVariableSolverData(fg, :x0)
@test all(isapprox.(mean(vnd.val), ProductRepr([0.0,0.0], [1.0 0.0; 0.0 1.0]), atol=0.1).parts)
@test all(is_point.(Ref(M), vnd.val))

vnd = getVariableSolverData(fg, :x1)
@test all(isapprox.(mean(vnd.val), ProductRepr([1.0,2.0], [0.7071 -0.7071; 0.7071 0.7071]), atol=0.1).parts)
@test all(is_point.(Ref(M), vnd.val))

v1 = addVariable!(fg, :x2, SpecialEuclidean2)
mf = ManifoldFactor(SpecialEuclidean(2), MvNormal([1,2,pi/4], [0.01,0.01,0.01]))
f = addFactor!(fg, [:x1, :x2], mf)

#test new error from solvetree
@test_broken solveTree!(fg; smtasks, verbose=true, recordcliqs=ls(fg)) isa Tuple

end


## ======================================================================================
##
## ======================================================================================
struct ManiPose2Point2{T <: SamplableBelief} <: IIF.AbstractManifoldMinimize
    Z::T
end

function IIF.getSample(cf::CalcFactor{<:ManiPose2Point2}, N::Int=1)
    ret = [rand(cf.factor.Z) for _ in 1:N]
    (ret, )
end

DFG.getManifold(::ManiPose2Point2) = TranslationGroup(2)

# define the conditional probability constraint
function (cfo::CalcFactor{<:ManiPose2Point2})(measX, p, q)
  #
    M = SpecialEuclidean(2)
    q_SE = ProductRepr(q, identity(SpecialOrthogonal(2), p.parts[2]))

    X_se2 = log(M, identity(M, p), compose(M, inv(M, p), q_SE))
    X = X_se2.parts[1]
    # NOTE wrong for what we want X̂ = log(M, p, q_SE)
    return measX - X 
end


@testset "Test SpecialEuclidean(2)" begin

Base.convert(::Type{<:Tuple}, M::TranslationGroup{Tuple{2},ℝ}) = (:Euclid, :Euclid)
Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{TranslationGroup{Tuple{2},ℝ}})  = (:Euclid, :Euclid)

@defVariable Point2 TranslationGroup(2) [0.0, 0.0]

##
fg = initfg()

v0 = addVariable!(fg, :x0, SpecialEuclidean2)

mp = ManifoldPrior(SpecialEuclidean(2), ProductRepr(@MVector([0.0,0.0]), @MMatrix([1.0 0.0; 0.0 1.0])), MvNormal([0.01, 0.01, 0.01]))
p = addFactor!(fg, [:x0], mp)

##
v1 = addVariable!(fg, :x1, Point2)
mf = ManiPose2Point2(MvNormal([1,2], [0.01,0.01]))
f = addFactor!(fg, [:x0, :x1], mf)


doautoinit!(fg, :x1)

vnd = getVariableSolverData(fg, :x1)
@test all(isapprox.(mean(vnd.val), [1.0,2.0], atol=0.1))

# ##
smtasks = Task[]
solveTree!(fg; smtasks, verbose=true, recordcliqs=ls(fg))
# # hists = fetchCliqHistoryAll!(smtasks);

vnd = getVariableSolverData(fg, :x0)
@test all(isapprox.(mean(vnd.val), ProductRepr([0.0,0.0], [1.0 0.0; 0.0 1.0]), atol=0.1).parts)

vnd = getVariableSolverData(fg, :x1)
@test all(isapprox.(mean(vnd.val), [1.0,2.0], atol=0.1))

end