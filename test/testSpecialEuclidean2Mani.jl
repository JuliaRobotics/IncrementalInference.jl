using DistributedFactorGraphs
using IncrementalInference
using Interpolations
using Manifolds
using StaticArrays
using Test
import IncrementalInference: HeatmapDensityRegular

## define new local variable types for testing

@defVariable Point2 TranslationGroup(2) [0.0, 0.0]

# @defVariable SpecialEuclidean2 SpecialEuclidean(2) ProductRepr(@MVector([0.0,0.0]), @MMatrix([1.0 0.0; 0.0 1.0]))
@defVariable SpecialEuclidean2 SpecialEuclidean(2) ProductRepr([0.0,0.0], [1.0 0.0; 0.0 1.0])

##

@testset "Test SpecialEuclidean(2)" begin

##

Base.convert(::Type{<:Tuple}, M::SpecialEuclidean{2}) = (:Euclid, :Euclid, :Circular)
Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{SpecialEuclidean{2}})  = (:Euclid, :Euclid, :Circular)


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


##

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
solveTree!(fg; smtasks, verbose=true) #, recordcliqs=ls(fg))
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

##

#test new error from solvetree
# smtasks = Task[]
@test solveTree!(fg; verbose=true) isa Tuple


## test partial prior issue

fg = initfg()

v0 = addVariable!(fg, :x0, SpecialEuclidean2)
mp = PartialPrior(MvNormal([0.01, 0.01]), (1,2))

p = addFactor!(fg, [:x0], mp, graphinit=false)

##

pbel_ = approxConvBelief(fg, :x0f1, :x0)

@test isPartial(pbel_)

@test pbel_._partial == [1;2]
@test length(pbel_.infoPerCoord) == 3

##
end

@testset "Test Pose2 like hex as SpecialEuclidean2" begin
##

M = getManifold(SpecialEuclidean2)
fg = initfg()
v0 = addVariable!(fg, :x0, SpecialEuclidean2)

mp = ManifoldPrior(SpecialEuclidean(2), ProductRepr(@MVector([10.0,10.0]), @MMatrix([-1.0 0.0; 0.0 -1.0])), MvNormal([0.1, 0.1, 0.01]))
p = addFactor!(fg, [:x0], mp)

##

for i in 0:5
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, SpecialEuclidean2)
    mf = ManifoldFactor(SpecialEuclidean(2), MvNormal([10.0,0,pi/3], [0.5,0.5,0.05]))
    f = addFactor!(fg, [psym;nsym], mf)
end


addVariable!(fg, :l1, SpecialEuclidean2, tags=[:LANDMARK;])
mf = ManifoldFactor(SpecialEuclidean(2), MvNormal([10.0,0,0], [0.1,0.1,0.01]))
addFactor!(fg, [:x0; :l1], mf)

mf = ManifoldFactor(SpecialEuclidean(2), MvNormal([10.0,0,0], [0.1,0.1,0.01]))
addFactor!(fg, [:x6; :l1], mf)

vnd = getVariableSolverData(fg, :x0)
@test isapprox(M, mean(M, vnd.val), ProductRepr([10.0,10.0], [-1.0 0.0; 0.0 -1.0]), atol=0.2)

vnd = getVariableSolverData(fg, :x1)
@test isapprox(M, mean(M, vnd.val), ProductRepr([0.0,10.0], [-0.5 0.866; -0.866 -0.5]), atol=0.4)

vnd = getVariableSolverData(fg, :x6)
@test isapprox(M, mean(M, vnd.val), ProductRepr([10.0,10.0], [-1.0 0.0; 0.0 -1.0]), atol=0.5)

##

smtasks = Task[]
solveTree!(fg; smtasks);

##
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
    q_SE = ProductRepr(q, identity_element(SpecialOrthogonal(2), p.parts[2]))

    X_se2 = log(M, identity_element(M, p), compose(M, inv(M, p), q_SE))
    X = X_se2.parts[1]
    # NOTE wrong for what we want X̂ = log(M, p, q_SE)
    return measX - X 
end


@testset "Test SpecialEuclidean(2)" begin

Base.convert(::Type{<:Tuple}, M::TranslationGroup{Tuple{2},ℝ}) = (:Euclid, :Euclid)
Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{TranslationGroup{Tuple{2},ℝ}})  = (:Euclid, :Euclid)

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


@testset "test propagateBelief w HeatmapSampler and init for PartialPriorPassThrough" begin
##

fg = initfg()

v0 = addVariable!(fg, :x0, SpecialEuclidean2)

img_ = rand(10,10).+5.0
x_,y_ = ([-9:2.0:9;],[-9:2.0:9;])

hmd = HeatmapDensityRegular(img_, (x_,y_), 5.5, 0.1, N=1000)
pthru = PartialPriorPassThrough(hmd, (1,2))

# test without nullhyp
f0 = addFactor!(fg, [:x0], pthru, graphinit=false)

## test the inference functions

bel, infd = propagateBelief(fg, v0, [f0])

@test isPartial(bel)
@test length(getPoints(bel)) == 1000


## repeat test with nullhypo

fg = initfg()

v0 = addVariable!(fg, :x0, SpecialEuclidean2)
# test with nullhypo
f0 = addFactor!(fg, [:x0], pthru, graphinit=false, nullhypo=0.2)

## test the inference functions

bel, infd = propagateBelief(fg, v0, [f0])
@test isPartial(bel)

## 

doautoinit!(fg, :x0)

@test length(getPoints(getBelief(fg, :x0))) == 1000
@info "PassThrough transfers the full point count to the graph, unless a product is calculated during the propagateBelief step."

##

solveGraph!(fg);

@test 1000 == length(getPoints(fg, :x0))

@warn "must still check if bandwidths are recalculated on many points (not necessary), or lifted from this case single prior"

##

mp = ManifoldPrior(SpecialEuclidean(2), ProductRepr(@MVector([0.0,0.0]), @MMatrix([1.0 0.0; 0.0 1.0])), MvNormal([0.01, 0.01, 0.01]))
f1 = addFactor!(fg, [:x0], mp, graphinit=false)

@test length(ls(fg, :x0)) == 2

##

prp, infd = propagateBelief(fg, v0, [f0;f1])

@test length(getPoints(prp)) == getSolverParams(fg).N

## check that solve corrects the point count on graph variable

@test 1000 == length(getPoints(fg, :x0))

solveGraph!(fg);

# this number should drop down to usual, 100 at time of writing
@test getSolverParams(fg).N == length(getPoints(fg, :x0))


##
end



#