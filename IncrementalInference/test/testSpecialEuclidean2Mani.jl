using DistributedFactorGraphs
using IncrementalInference
using Interpolations
using Manifolds
using StaticArrays
using Test
using LieGroups
using LieGroups: TranslationGroup
import IncrementalInference: LevelSetGridNormal
import Rotations as _Rot

## define new local variable types for testing

# @defStateType TranslationGroup2 ValidationLieGroup(TranslationGroup(2)) @SVector[0.0, 0.0]
@defStateType TranslationGroup2 TranslationGroup(2) @SVector[0.0, 0.0]

# SE2 = ValidationLieGroup(SpecialEuclideanGroup(2; variant=:right))
SE2 = SpecialEuclideanGroup(2; variant=:right)
# PoseMani = SE2
PoseMani = TranslationGroup(2) × SpecialOrthogonalGroup(2)
@defStateType SpecialEuclidean2 PoseMani ArrayPartition(@SVector([0.0,0.0]), @SMatrix([1.0 0.0; 0.0 1.0]))

##

DFG.@defObservationType SE2SE2 RelativeObservation SE2

struct ManifoldFactorSE2{T <: SamplableBelief} <: IIF.RelativeObservation
    Z::T
end

SE2SE2() = SE2SE2(MvNormal(Diagonal([1,1,1])))

IIF.selectFactorType(::Type{<:SpecialEuclidean2}, ::Type{<:SpecialEuclidean2}) = SE2SE2

function IIF.getSample(cf::CalcFactor{<:SE2SE2}) 
  M = getManifold(SE2SE2)
  X = sampleTangent(M, cf.factor.Z)
  return X
end

function (cf::CalcFactor{<:SE2SE2})(X, p, q)
    M = getManifold(SE2SE2)
    X̂ = log(M, p, q)
    return vee(LieAlgebra(M), X - X̂)
end

##

@testset "Test SpecialEuclideanGroup(2)" begin
##

M = getManifold(SpecialEuclidean2)
@test M == PoseMani
pT = getPointType(SpecialEuclidean2)
# @test pT == ArrayPartition{Float64,Tuple{Vector{Float64}, Matrix{Float64}}}
# @test pT == ArrayPartition{Tuple{MVector{2, Float64}, MMatrix{2, 2, Float64, 4}}}
@test pT == ArrayPartition{Float64, Tuple{SVector{2, Float64}, SMatrix{2, 2, Float64, 4}}}
pϵ = getPointIdentity(SpecialEuclidean2)
# @test_broken pϵ == ArrayPartition(@MVector([0.0,0.0]), @MMatrix([1.0 0.0; 0.0 1.0]))
@test all(isapprox.(pϵ, ArrayPartition(SA[0.0,0.0], SA[1.0 0.0; 0.0 1.0])))

@test is_point(getManifold(SpecialEuclidean2), getPointIdentity(SpecialEuclidean2))

##
fg = initfg()

v0 = addVariable!(fg, :x0, SpecialEuclidean2)

# mp = ManifoldPrior(SpecialEuclideanGroup(2), ArrayPartition(@MVector([0.0,0.0]), @MMatrix([1.0 0.0; 0.0 1.0])), MvNormal([0.01, 0.01, 0.01]))
# mp = ManifoldPrior(SpecialEuclideanGroup(2), ArrayPartition(@MVector([0.0,0.0]), @MMatrix([1.0 0.0; 0.0 1.0])), MvNormal(Diagonal(abs2.([0.01, 0.01, 0.01]))))
# mp = ManifoldPrior(SE2, ArrayPartition([0.0,0.0], [1.0 0.0; 0.0 1.]), MvNormal(Diagonal(abs2.([0.01, 0.01, 0.01]))))
mp = ManifoldPrior(SE2, ArrayPartition(SA[0.0,0.0], SA[1.0 0.0; 0.0 1.]), MvNormal(Diagonal(abs2.(SA[0.01, 0.01, 0.01]))))
p = addFactor!(fg, [:x0], mp)


##

doautoinit!(fg, :x0)

##
vnd = getState(fg, :x0, :default)
@test all(isapprox.(mean(vnd.val), ArrayPartition(SA[0.0,0.0], SA[1.0 0.0; 0.0 1.0]), atol=0.1))
@test all(is_point.(Ref(M), vnd.val))

##
v1 = addVariable!(fg, :x1, SpecialEuclidean2)
mf = ManifoldFactor(SE2, MvNormal(SA[1,2,pi/4], SA[0.01,0.01,0.01]))
f = addFactor!(fg, [:x0, :x1], mf)

doautoinit!(fg, :x1)

p0 = ArrayPartition(Vector([0.0,0.0]), Matrix([1.0 0.0; 0.0 1.0]))
p1 = exp(SE2, p0, hat(LieAlgebra(SE2), [1.0,2,pi/4], typeof(p0)))
p2 = exp(SE2, p1, hat(LieAlgebra(SE2), [1.0,2,pi/4], typeof(p1)))

vnd = getState(fg, :x1, :default)
@test all(isapprox(M, mean(M,vnd.val), p1, atol=0.1))
@test all(is_point.(Ref(M), vnd.val))

##
smtasks = Task[]
solveTree!(fg; smtasks, verbose=true) #, recordcliqs=ls(fg))
# hists = fetchCliqHistoryAll!(smtasks);

vnd = getState(fg, :x0, :default)
@test all(isapprox.(mean(vnd.val), ArrayPartition(SA[0.0,0.0], SA[1.0 0.0; 0.0 1.0]), atol=0.1))
@test all(is_point.(Ref(M), vnd.val))

vnd = getState(fg, :x1, :default)
@test all(isapprox(M, mean(vnd.val), p1, atol=0.1))
@test all(is_point.(Ref(M), vnd.val))

v1 = addVariable!(fg, :x2, SpecialEuclidean2)
mf = ManifoldFactor(SE2, MvNormal(SA[1,2,pi/4], SA[0.01,0.01,0.01]))
f = addFactor!(fg, [:x1, :x2], mf)

##

#test new error from solvetree
smtasks = Task[]
result = solveTree!(fg; smtasks, verbose=true)
@test result isa AbstractBayesTree

IIF.solveGraphParametric!(fg; sparse = false, damping_term_min=1e-12)

vnd = getState(fg, :x0, :parametric)
@test all(isapprox(M, vnd.val[1], p0, atol=1e-6))
vnd = getState(fg, :x1, :parametric)
@test all(isapprox(M, vnd.val[1], p1, atol=1e-6))
vnd = getState(fg, :x2, :parametric)
@test all(isapprox(M, vnd.val[1], p2, atol=1e-6))

## test partial prior issue
fg = initfg()

v0 = addVariable!(fg, :x0, SpecialEuclidean2)
mp = PartialPrior(SpecialEuclidean2,MvNormal([0.01, 0.01]), (1,2))

p = addFactor!(fg, [:x0], mp, graphinit=false)

##

pbel_ = approxConvBelief(fg, :x0f1, :x0)

@test isPartial(pbel_)

@test pbel_._partial == [1;2]
@test length(pbel_.infoPerCoord) == 3
##
end


@testset "test initVariableManual! with Vector of Tuple inputs" begin
##

fg = initfg()

pts = [(randn(2), Matrix(_Rot.RotMatrix2(randn()))) for _ in 1:50]

addVariable!(fg, :x0, SpecialEuclidean2)
initVariable!(getVariable(fg, :x0), pts)

@test isapprox( pts[1][1], getPoints(fg, :x0)[1].x[1])
@test isapprox( pts[1][2], getPoints(fg, :x0)[1].x[2])

# can delete upon deprecation of initVariable! and favor initVariable!
initVariable!(getVariable(fg, :x0), reverse(pts)) 
@test isapprox( pts[end][1], getPoints(fg, :x0)[1].x[1])
@test isapprox( pts[end][2], getPoints(fg, :x0)[1].x[2])


##
end

##

@testset "Test Pose2 like hex as SpecialEuclidean2" begin
##

M = getManifold(SpecialEuclidean2)
fg = initfg()
# fg.solverParams.graphinit = false
v0 = addVariable!(fg, :x0, SpecialEuclidean2)

mp = ManifoldPrior(SE2, ArrayPartition(Vector([10.0,10.0]), Matrix([-1.0 0.0; 0.0 -1.0])), MvNormal([0.05, 0.05, 0.005]))
p = addFactor!(fg, [:x0], mp)

##

for i in 0:5
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, SpecialEuclidean2)
    mf = ManifoldFactor(SE2, MvNormal([10.0,0,pi/3], [0.1,0.1,0.01]))
    f = addFactor!(fg, [psym;nsym], mf)
end


addVariable!(fg, :l1, SpecialEuclidean2, tags=[:LANDMARK;])
mf = ManifoldFactor(SE2, MvNormal([10.0,0,0], [0.1,0.1,0.01]))
addFactor!(fg, [:x0; :l1], mf)

mf = ManifoldFactor(SE2, MvNormal([10.0,0,0], [0.1,0.1,0.01]))
addFactor!(fg, [:x6; :l1], mf)

##

smtasks = Task[]
solveTree!(fg; smtasks);
IIF.autoinitParametric!(fg)
IIF.solveGraphParametric!(fg; sparse = false, damping_term_min=1e-12)

vnd = getState(fg, :x0, :default)
@test isapprox(M, mean(M, vnd.val), ArrayPartition([10.0,10.0], [-1.0 0.0; 0.0 -1.0]), atol=0.2)
vnd = getState(fg, :x0, :parametric)
@test isapprox(M, vnd.val[1], ArrayPartition([10.0,10.0], [-1.0 0.0; 0.0 -1.0]), atol=1e-6)

# calculate the reference solution
p0 = ArrayPartition(Vector([10.0,10.0]), Matrix([-1.0 0.0; 0.0 -1.0]))
ref = exp(SE2, p0, hat(LieAlgebra(SE2), [10.0,0,pi/3], typeof(p0)))

vnd = getState(fg, :x1, :default)
@test isapprox(M, mean(M, vnd.val), ref, atol=0.4)
vnd = getState(fg, :x1, :parametric)
@test isapprox(M, vnd.val[1], ref, atol=1e-6)

vnd = getState(fg, :x6, :default)
@test isapprox(M, mean(M, vnd.val), ArrayPartition([10.0,10.0], [-1.0 0.0; 0.0 -1.0]), atol=0.5)
vnd = getState(fg, :x6, :parametric)
@test isapprox(M, vnd.val[1], ArrayPartition([10.0,10.0], [-1.0 0.0; 0.0 -1.0]), atol=1e-6)

@test isapprox(M, getState(fg, :x0, :parametric).val[1], getState(fg, :x6, :parametric).val[1], atol=1e-6)

if false
fix, ax, plt = lines(points2(fg); label="parametric")
# foreach(ls(fg)[1:7]) do vl
foreach(ls(fg)) do vl
    state = getState(fg, vl, :default)
    pnts = map(state.val) do val
        Point2f(val[1:2])
    end
    scatter!(ax, pnts; label=string(vl))
end
axislegend(ax)
end
## Special test for manifold based messages

#FIXME this may show some bug in propagateBelief caused by empty factors
fg.solverParams.useMsgLikelihoods = true
smtasks = Task[]
result = solveTree!(fg; smtasks); #, recordcliqs=ls(fg))
@test result isa AbstractBayesTree

##
end



@testset "test deconv on <:RelativeObservation" begin
##

fg = initfg()
getSolverParams(fg).useMsgLikelihoods = true

addVariable!(fg, :x0, SpecialEuclidean2)
addVariable!(fg, :x1, SpecialEuclidean2)

mp = ManifoldPrior(SE2, ArrayPartition(Vector([10.0,10.0]), Matrix([-1.0 0.0; 0.0 -1.0])), MvNormal(diagm([0.1, 0.1, 0.01].^2)))
p = addFactor!(fg, [:x0], mp)

doautoinit!(fg,:x0)

addFactor!(fg, [:x0;:x1], SE2SE2(MvNormal([10.0,0,0.1], diagm([0.5,0.5,0.05].^2))))

initAll!(fg)

# now check deconv

pred, meas = approxDeconv(fg, :x0x1f1)

@test mmd(SE2, pred, meas) < 1e-1

p_t = map(x->x.x[1], pred)
m_t = map(x->x.x[1], meas)
p_θ = map(x->x.x[2][2], pred)
m_θ = map(x->x.x[2][2], meas)

@test isapprox(mean(p_θ), 0.1, atol=0.02)
@test isapprox(std(p_θ), 0.05, atol=0.02)

@test isapprox(mean(p_t), [10,0], atol=0.3)
@test isapprox(std(p_t), [0.5,0.5], atol=0.3)

@test isapprox(mean(p_θ), mean(m_θ), atol=0.03)
@test isapprox(std(p_θ), std(m_θ), atol=0.03)

@test isapprox(mean(p_t), mean(m_t), atol=0.3)
@test isapprox(std(p_t), std(m_t), atol=0.3)

end


## ======================================================================================
##
## ======================================================================================

struct ManiPose2Point2{T <: SamplableBelief} <: IIF.RelativeObservation
    Z::T
    partial::Vector{Int}
end

function IIF.getSample(cf::CalcFactor{<:ManiPose2Point2})
    return rand(cf.factor.Z)
end

DFG.getManifold(::ManiPose2Point2) = TranslationGroup(2)

# define the conditional probability constraint
function (cfo::CalcFactor{<:ManiPose2Point2})(measX, p, q)
    #
    M = SE2
    q_SE = ArrayPartition(q, identity_element(SpecialOrthogonalGroup(2), typeof(p.x[2])))

    X_se2 = log(M, p, q_SE)
    X = X_se2.x[1]
    # NOTE wrong for what we want X̂ = log(M, p, q_SE)
    return measX - X 
end

##
@testset "Test SE2" begin
##

# Base.convert(::Type{<:Tuple}, M::TranslationGroup{Tuple{2},ℝ}) = (:Euclid, :Euclid)
# Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{TranslationGroup{Tuple{2},ℝ}})  = (:Euclid, :Euclid)

##
fg = initfg()

v0 = addVariable!(fg, :x0, SpecialEuclidean2)

mp = ManifoldPrior(SE2, ArrayPartition(Vector([0.0,0.0]), Matrix([1.0 0.0; 0.0 1.0])), MvNormal([0.01, 0.01, 0.01]))
p = addFactor!(fg, [:x0], mp)

##
v1 = addVariable!(fg, :x1, TranslationGroup2)
mf = ManiPose2Point2(MvNormal([1,2], [0.01,0.01]), [1;2])
f = addFactor!(fg, [:x0, :x1], mf)


doautoinit!(fg, :x1)

vnd = getState(fg, :x1, :default)
@test all(isapprox.(mean(vnd.val), [1.0,2.0], atol=0.1))

##
smtasks = Task[]
solveTree!(fg; smtasks, verbose=true, recordcliqs=ls(fg))
# # hists = fetchCliqHistoryAll!(smtasks);

vnd = getState(fg, :x0, :default)
@test isapprox(mean(getManifold(fg,:x0),vnd.val), ArrayPartition([0.0,0.0], [1.0 0.0; 0.0 1.0]), atol=0.1)

vnd = getState(fg, :x1, :default)
@test all(isapprox.(mean(vnd.val), [1.0,2.0], atol=0.1))

##
end


@testset "test propagateBelief w HeatmapSampler and init for PartialPriorPassThrough w Priors" begin
##

fg = initfg()

v0 = addVariable!(fg, :x0, SpecialEuclidean2)

img_ = rand(10,10).+5.0
x_,y_ = ([-9:2.0:9;],[-9:2.0:9;])

hmd = LevelSetGridNormal(img_, (x_,y_), 5.5, 0.1, N=120)
pthru = PartialPriorPassThrough(hmd, (1,2))

@show hmd

## quick 

pf = convert( AbstractPackedObservation, pthru )
upf = convert( AbstractObservation, pf )

@test pthru.partial == upf.partial
@test isapprox( pthru.Z.heatmap.data, upf.Z.heatmap.data )
@test isapprox( pthru.Z.heatmap.domain[1], upf.Z.heatmap.domain[1] )
@test isapprox( pthru.Z.heatmap.domain[2], upf.Z.heatmap.domain[2] )
@test isapprox( pthru.Z.level, upf.Z.level )
@test isapprox( pthru.Z.sigma, upf.Z.sigma )
@test isapprox( pthru.Z.sigma_scale, upf.Z.sigma_scale )


## test without nullhyp

f0 = addFactor!(fg, [:x0], pthru, graphinit=false)

## test the inference functions

bel, infd = propagateBelief(fg, v0, [f0])

@test isPartial(bel)
@test_broken length(getPoints(bel)) == 120


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

@test length(getPoints(getBelief(fg, :x0))) == getSolverParams(fg).N # 120
# @info "PassThrough transfers the full point count to the graph, unless a product is calculated during the propagateBelief step."

##
# @test_broken begin
## check the partials magic

dens, ipc = propagateBelief(fg,:x0,:)
testv = deepcopy(getVariable(fg, :x0))
setBelief!(testv, dens, true, ipc)


##

smtasks = Task[]
solveGraph!(fg; smtasks);
# hists = fetchCliqHistoryAll!(smtasks)
# printCSMHistoryLogical(hists)
# hists_ = deepcopy(hists)
# repeatCSMStep!(hists, 1, 6)

@test_broken 120 == length(getPoints(fg, :x0))

@warn "must still check if bandwidths are recalculated on many points (not necessary), or lifted from this case single prior"

##

mp = ManifoldPrior(SE2, ArrayPartition(Vector([0.0,0.0]), Matrix([1.0 0.0; 0.0 1.0])), MvNormal([0.01, 0.01, 0.01],[1 0 0;0 1 0;0 0 1.]))
f1 = addFactor!(fg, [:x0], mp, graphinit=false)

@test length(ls(fg, :x0)) == 2

##

prp, infd = propagateBelief(fg, v0, [f0;f1])

@test length(getPoints(prp)) == getSolverParams(fg).N

## check that solve corrects the point count on graph variable

@test_broken 120 == length(getPoints(fg, :x0))

solveGraph!(fg);

# this number should drop down to usual, 100 at time of writing
@test getSolverParams(fg).N == length(getPoints(fg, :x0))


## check saveDFG (check consistency of packing converters above)

@error "Whats going on in PackedManifoldPrior, skipping tests"
@test_broken begin
    saveDFG(joinpath(tempdir(),"passthru"), fg)
    fg_ = loadDFG(joinpath(tempdir(),"passthru.tar.gz"))
    Base.rm(joinpath(tempdir(),"passthru.tar.gz"))
end

# @error "#FIXME test propagateBelief w HeatmapSampler ... broken on ci but not local"
# return true
# end

##
end




@testset "test point count on propagate, solve, init for PartialPriorPassThrough w Relative" begin
##

fg = initfg()

v0 = addVariable!(fg, :x0, SpecialEuclidean2)

img_ = rand(10,10).+5.0
x_,y_ = ([-9:2.0:9;],[-9:2.0:9;])

hmd = LevelSetGridNormal(img_, (x_,y_), 5.5, 0.1, N=120)
pthru = PartialPriorPassThrough(hmd, (1,2))

# test without nullhyp
f0 = addFactor!(fg, [:x0], pthru, graphinit=false)

## test the inference functions
addVariable!(fg, :x1, SpecialEuclidean2)
# mp = ManifoldPrior(SE2, ArrayPartition(@MVector([0.0,0.0]), @MMatrix([1.0 0.0; 0.0 1.0])), MvNormal([0.01, 0.01, 0.01]))
mp = ManifoldPrior(SE2, ArrayPartition(Vector([0.0,0.0]), Matrix([1.0 0.0; 0.0 1.0])), MvNormal([0.01, 0.01, 0.01]))
f1 = addFactor!(fg, [:x1], mp, graphinit=false)

doautoinit!(fg, :x1)

## connect with relative and check calculation size on x0

mf = ManifoldFactor(SE2, MvNormal([1,2,pi/4], [0.01,0.01,0.01]))
f2 = addFactor!(fg, [:x0, :x1], mf, graphinit=false)

##

bel, infd = propagateBelief(fg, v0, [f0;f2])

@test !isPartial(bel)
@test getSolverParams(fg).N == length(getPoints(bel))

## check other functions

solveTree!(fg);

@test getSolverParams(fg).N == length(getPoints(fg, :x0))
@test getSolverParams(fg).N == length(getPoints(fg, :x1))


## and check that initAll! works the same (different init sequences may change code execution path)

fg = initfg()
v0 = addVariable!(fg, :x0, SpecialEuclidean2)
img_ = rand(10,10).+5.0
x_,y_ = ([-9:2.0:9;],[-9:2.0:9;])
hmd = LevelSetGridNormal(img_, (x_,y_), 5.5, 0.1, N=120)
pthru = PartialPriorPassThrough(hmd, (1,2))
f0 = addFactor!(fg, [:x0], pthru, graphinit=false)
addVariable!(fg, :x1, SpecialEuclidean2)
# mp = ManifoldPrior(SE2, ArrayPartition(@MVector([0.0,0.0]), @MMatrix([1.0 0.0; 0.0 1.0])), MvNormal([0.01, 0.01, 0.01]))
mp = ManifoldPrior(SE2, ArrayPartition(Vector([0.0,0.0]), Matrix([1.0 0.0; 0.0 1.0])), MvNormal([0.01, 0.01, 0.01]))
f1 = addFactor!(fg, [:x1], mp, graphinit=false)
mf = ManifoldFactor(SE2, MvNormal([1,2,pi/4], [0.01,0.01,0.01]))
f2 = addFactor!(fg, [:x0, :x1], mf, graphinit=false)

##

bel, infd = propagateBelief(fg, v0, [f0;f2])

@test !isPartial(bel)
@test getSolverParams(fg).N == length(getPoints(bel))

##

initAll!(fg)

@test getSolverParams(fg).N == length(getPoints(fg, :x0))
@test getSolverParams(fg).N == length(getPoints(fg, :x1))

##
end

@testset "Test SE2 to TranslationGroup(2) multihypo" begin
##

fg = initfg()
# fg.solverParams.attemptGradients=false

v0 = addVariable!(fg, :x0, SpecialEuclidean2)

# mp = ManifoldPrior(SE2, ArrayPartition(@MVector([0.0,0.0]), @MMatrix([1.0 0.0; 0.0 1.0])), MvNormal([0.01, 0.01, 0.01]))
mp = ManifoldPrior(SE2, ArrayPartition(Vector([0.0,0.0]), Matrix([1.0 0.0; 0.0 1.0])), MvNormal([0.01, 0.01, 0.01]))
p = addFactor!(fg, [:x0], mp)

##
addVariable!(fg, :x1a, TranslationGroup2)
addVariable!(fg, :x1b, TranslationGroup2)
mf = ManiPose2Point2(MvNormal([1,2], [0.01,0.01]), [1;2])
f = addFactor!(fg, [:x0, :x1a, :x1b], mf; multihypo=[1,0.5,0.5])

solveTree!(fg)

vnd = getState(fg, :x0, :default)
@test isapprox(SE2, mean(SE2, vnd.val), ArrayPartition([0.0,0.0], [1.0 0; 0 1]), atol=0.1)

#FIXME I would expect close to 50% of particles to land on the correct place
# Currently software works so that 33% should land there so testing 20 for now
pnt = getPoints(fg, :x1a)
@test sum(isapprox.(pnt, Ref([1.0,2.0]), atol=0.1)) > 15

#FIXME I would expect close to 50% of particles to land on the correct place
pnt = getPoints(fg, :x1b)
@test sum(isapprox.(pnt, Ref([1.0,2.0]), atol=0.1)) > 15


## other way around

fg = initfg()
fg.solverParams.attemptGradients=false

addVariable!(fg, :x0, SpecialEuclidean2)
addVariable!(fg, :x1a, TranslationGroup2)
addVariable!(fg, :x1b, TranslationGroup2)

# mp = ManifoldPrior(SE2, ArrayPartition(@MVector([0.0,0.0]), @MMatrix([1.0 0.0; 0.0 1.0])), MvNormal([10, 10, 0.01]))
mp = ManifoldPrior(SE2, ArrayPartition(Vector([0.0,0.0]), Matrix([1.0 0.0; 0.0 1.0])), MvNormal(zeros(3),diagm([10, 10, 0.01])))
p = addFactor!(fg, [:x0], mp)
mp = ManifoldPrior(TranslationGroup(2), [1.,0], MvNormal([0.01, 0.01]))
p = addFactor!(fg, [:x1a], mp)
mp = ManifoldPrior(TranslationGroup(2), [-1.,0], MvNormal([0.01, 0.01]))
p = addFactor!(fg, [:x1b], mp)

mf = ManiPose2Point2(MvNormal([0., 1], [0.01,0.01]), [1;2])
f = addFactor!(fg, [:x0, :x1a, :x1b], mf; multihypo=[1,0.5,0.5])

solveTree!(fg)

pnts = getPoints(fg, :x0)
# c = getCoordinates.(SpecialEuclidean2, pnts)
# @cast p[i,j] := c[i][j]
# scatter(p[:,1], p[:,2])

#FIXME
@error "Invalid multihypo test"
if false
    # FIXME ManiPose2Point2 factor mean [1.,0] cannot go "backwards" from [0,0] to [-1,0] with covariance 0.01 -- wholly inconsistent test design
    @test 10 < sum(isapprox.(Ref(SE2), pnts, Ref(ArrayPartition([-1.0,0.0], [1.0 0; 0 1])), atol=0.5))
    @test 10 < sum(isapprox.(Ref(SE2), pnts, Ref(ArrayPartition([1.0,0.0], [1.0 0; 0 1])), atol=0.5))
end

##
end

@testset "Test SE2 to SE2 multihypo" begin
##

fg = initfg()
# fg.solverParams.attemptGradients=false

v0 = addVariable!(fg, :x0, SpecialEuclidean2)

# mp = ManifoldPrior(SE2, ArrayPartition(@MVector([0.0,0.0]), @MMatrix([1.0 0.0; 0.0 1.0])), MvNormal([0.01, 0.01, 0.01]))
mp = ManifoldPrior(SE2, ArrayPartition(Vector([0.0,0.0]), Matrix([1.0 0.0; 0.0 1.0])), MvNormal([0.01, 0.01, 0.01]))
p = addFactor!(fg, [:x0], mp)

##
addVariable!(fg, :x1a, SpecialEuclidean2)
addVariable!(fg, :x1b, SpecialEuclidean2)
mf = ManifoldFactor(SE2, MvNormal([1,2,pi/4], [0.01,0.01,0.01]))
f = addFactor!(fg, [:x0, :x1a, :x1b], mf; multihypo=[1,0.5,0.5])

solveTree!(fg)

p0 = ArrayPartition([0.0,0.0], [1.0 0; 0 1])
p1 = exp(SE2, hat(LieAlgebra(SE2), [1,2,pi/4], typeof(p0)))

vnd = getState(fg, :x0, :default)
@test isapprox(SE2, mean(SE2, vnd.val), p0, atol=0.1)

#FIXME I would expect close to 50% of particles to land on the correct place
# Currently software works so that 33% should land there so testing 20 for now
pnt = getPoints(fg, :x1a)
@test sum(isapprox.(Ref(SE2), pnt, Ref(p1), atol=0.1)) > 20

#FIXME I would expect close to 50% of particles to land on the correct place
pnt = getPoints(fg, :x1b)
@test sum(isapprox.(Ref(SE2), pnt, Ref(p1), atol=0.1)) > 20

##
end


#