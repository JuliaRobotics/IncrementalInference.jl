# partial constraint development
using Test
using IncrementalInference
# going to introduce two new constraint mutable structs
using Statistics
using TensorCast
using Manifolds

import IncrementalInference: getManifold, getSample

##

mutable struct DevelopPartial{P <: Tuple} <: AbstractPrior
  x::Distribution
  partial::P 
end
getSample(cf::CalcFactor{<:DevelopPartial}) = rand(cf.factor.x, 1)

#
mutable struct DevelopDim2 <: AbstractPrior
  x::Distribution
end
getSample(cf::CalcFactor{<:DevelopDim2}) = rand(cf.factor.x, 1)


mutable struct DevelopPartialPairwise <: AbstractRelativeMinimize
  x::Distribution
  partial::Tuple
  DevelopPartialPairwise(x::Distribution) = new(x, (2,))
end
getManifold(::IIF.InstanceType{DevelopPartialPairwise}) = TranslationGroup(1)

getSample(cf::CalcFactor{<:DevelopPartialPairwise}) = rand(cf.factor.x, 1)

function (dp::CalcFactor{<:DevelopPartialPairwise})(meas,
                                                    x1,
                                                    x2  )
  #
  return meas[1] - (x2[2]-x1[2])
end


## prep

N=100

pr = DevelopDim2(MvNormal([0.0;0.0], diagm([0.01;0.01])))
dp = DevelopPartial(Normal(2.0, 1.0),(1,))

## build graph

fg = initfg()

v1 = addVariable!(fg,:x1,ContinuousEuclid{2}(),N=N)
f1  = addFactor!(fg,[:x1],pr)
f2  = addFactor!(fg,[:x1],dp, graphinit=false)

doautoinit!(fg, :x1)

##

@testset "test evaluation of full constraint prior" begin


pts_, _ = evalFactor(fg, f1, v1.label, N=N)
@cast pts[i,j] := pts_[j][i]
@test size(pts,1) == 2
@test size(pts,2) == N
@test norm(Statistics.mean(pts,dims=2)[1] .- [0.0]) < 0.3


end

##

@testset "test evaluation of partial constraint prior" begin
##

memcheck_ = getVal(v1)
@cast memcheck[i,j] := memcheck_[j][i]

X1pts_ = getVal(v1)
@cast X1pts[i,j] := X1pts_[j][i]
pts_ = approxConv(fg, getLabel(f2), :x1, N=N)
@cast pts[i,j] := pts_[j][i]

@test size(pts, 1) == 2
@test size(pts, 2) == N
@test norm(Statistics.mean(pts,dims=2)[1] .- [2.0]) < 0.75
# ensure the correct response from
@test norm(X1pts[1,:] - pts[1,:]) > 2.0
@test norm(X1pts[2,:] - pts[2,:]) < 1e-10
@test norm(X1pts - memcheck) < 1e-10

##
end


@testset "check that partials are received through convolutions of prior" begin
##

# check that a partial belief is obtained
X1_ = approxConvBelief(fg, :x1f2, :x1)

@test isPartial(X1_)

##
end


@testset "test solving of factor graph" begin
##

getSolverParams(fg).N = N
tree = solveTree!(fg)
pts_ = getVal(fg, :x1)
@cast pts[i,j] := pts_[j][i]

@test norm(Statistics.mean(pts,dims=2)[1] .- [0.0]) < 0.4
@test norm(Statistics.mean(pts,dims=2)[2] .- [0.0]) < 0.4

##
end
# plotKDE(getBelief(fg, :x1),levels=3)


## partial relative gradient and graph

v2 = addVariable!(fg,:x2,ContinuousEuclid{2},N=N)

dpp = DevelopPartialPairwise(Normal(10.0, 1.0))
f3  = addFactor!(fg,[:x1;:x2],dpp)

dp2 = DevelopPartial( Normal(-20.0, 1.0), (1,) )
f4  = addFactor!(fg,[:x2;], dp2, graphinit=false)

doautoinit!(fg, :x2)

##

@testset "test partial info per coord through relative convolution (conditional)" begin
##

one_meas = [10.0;]
pts = ([0;0.0], [0;10.0])
gradients = FactorGradientsCached!(dpp, (ContinuousEuclid{2}, ContinuousEuclid{2}), one_meas, pts);

##

# check that the gradients can be calculated
J = gradients(one_meas, pts...)

@test size(J) == (4,4)
@test norm(J - [0 0 0 0; 0 0 0 1; 0 0 0 0; 0 1 0 0] ) < 1e-4

## check perturbation logic

prtb = calcPerturbationFromVariable(gradients, [1=>[1;1]])

# self variation is taken as 0 at this time
@test isapprox( prtb[1], [0;0] )
# variable 1 influences 2 only through partial dimension 2 (as per DevelopPartialPairwise)
@test isapprox( prtb[2], [0;1] )

##  test evaluation through the convolution operation withing a factor graph

# add relative IPC calculation inside evalFactor
bel = approxConvBelief(fg, getLabel(f3), :x2)
@test isPartial(bel)

##
end

##

@testset "test evaluation of multiple simultaneous partial constraints" begin
global fg

##

initAll!(fg)
valx2_ = getVal(fg, :x2)
@cast valx2[i,j] := valx2_[j][i]
pts_ = approxConv(fg, :x1x2f1, :x2, N=N)
@cast pts[i,j] := pts_[j][i]
@test size(pts,1) == 2
@test norm(Statistics.mean(pts,dims=2)[2] .- [10.0]) < 3.0
@test norm(valx2[1,:] - pts[1,:]) < 1e-5

pts_ = approxConv(fg, :x2f1, :x2, N=N)
@cast pts[i,j] := pts_[j][i]
@test size(pts,1) == 2
@test norm(Statistics.mean(pts,dims=2)[1] .- [-20.0]) < 0.75
@test (Statistics.std(pts,dims=2)[1] .- 1.0) < 0.4

##

end

##

# keep previous values to ensure funciton evaluation is modifying correct data fields

@warn "restore calcProposalBelief as testset!"
# @testset "test calcProposalBelief..." begin
# global v2, fg, f3, f4, N


thefac = getFactor(fg, :x1x2f1)

X2lpts_ = getVal(getVariable(fg, :x2))
@cast X2lpts[i,j] := X2lpts_[j][i]
keepaside, = (calcProposalBelief(fg, thefac, :x2, N=N),)
@test Ndim(keepaside) == 2
lpts_ = getPoints(keepaside, false)
@cast lpts[i,j] := lpts_[j][i]
@test length(lpts_) == N

@show X2lpts[2,95:100]
@show lpts[2,95:100]
@show getPoints(keepaside)

##

# DevelopPartialPairwise must only modify the second dimension of proposal distribution on X2
@test norm(X2lpts[1,:] - lpts[1,:]) < 1e-10
# @test norm(X2lpts[2,:] - lpts[2,:]) > 1e-10 # 10*N*0.5 # why so big?
memcheck_ = getVal(v2)
@cast memcheck[i,j] := memcheck_[j][i]
@test norm(X2lpts - memcheck) < 1e-10


X2lpts_ = getVal(v2)
@cast X2lpts[i,j] := X2lpts_[j][i]
p4 = calcProposalBelief(fg, f4, v2.label, N=N)
@test Ndim(p4) == 2
lpts_ = getPoints(keepaside, false)
@cast lpts[i,j] := lpts_[j][i]
@test length(lpts_) == N

# DevelopPartialPairwise must only modify the second dimension of proposal distribution on X2
@test norm(X2lpts[1,:] - lpts[1,:]) < 1e-10
@test norm(X2lpts[2,:] - lpts[2,:]) > 1e-10 # 10*N*0.5 # why so big?
memcheck_ = getVal(v2)
@cast memcheck[i,j] := memcheck_[j][i]
@test norm(X2lpts - memcheck) < 1e-10

# end

##


@testset "test belief prediction with partials..." begin

##

global v2, fg

# partial prior
X2pts_ = getVal(v2)
@cast X2pts[i,j] := X2pts_[j][i]
# NOTE, SUPER IMPORTANT, predictbelief returns full dimension points (even if only partials are sent in for proposals)
val_, = predictbelief(fg, v2, [f4], N=N)
@cast val[i,j] := val_[j][i]
@show X2pts_[1]';
@show val_[1]';
@test norm(X2pts[2,:] - val[2,:]) < 1e-10
@test 0.0 < norm(X2pts[1,:] - val[1,:])
@test norm(Statistics.mean(val[1,:]) .+ 20.0) < 0.75


# partial pairwise
X2pts_ = getVal(v2)
@cast X2pts[i,j] := X2pts_[j][i]
val_, = predictbelief(fg, v2, [f3], N=N)
@cast val[i,j] := val_[j][i]
@test norm(X2pts[1,:] - val[1,:]) < 1e-10
@test 0.0 < norm(X2pts[2,:] - val[2,:])
val2_ = getVal(v1)
@cast val2[i,j] := val2_[j][i]
@test abs(Statistics.mean(val[2,:] - val2[2,:]) .- 10.0) < 0.75

##

# combination of partials
val_, = predictbelief(fg, v2, [f3;f4], N=N)
@cast val[i,j] := val_[j][i]
# plotKDE(kde!(val),levels=3)
@test norm(Statistics.mean(val,dims=2)[1] .- [-20.0]) < 1
@error "restore inconsistent test result (not always broken)"
if false
  @test_broken norm(Statistics.mean(val,dims=2)[2] .- [10.0]) < 0.01
end
@test (Statistics.std(val,dims=2)[1] .- 1.0) < 3.0
@test (Statistics.std(val,dims=2)[2] .- 1.0) < 3.0

##

getSolverParams(fg).N = N
tree = solveTree!(fg)


pts_ = getVal(fg, :x1)
@cast pts[i,j] := pts_[j][i]
@test norm(Statistics.mean(pts,dims=2)[1] .- [0.0]) < 0.5
@test norm(Statistics.mean(pts,dims=2)[2] .- [0.0]) < 0.5

pts_ = getVal(fg, :x2)

ppe = getPPE(fg, :x2).mean

X2 = getBelief(fg, :x2)

# check mean is close
@test_broken isapprox(mean(X2), [-20;10], atol=0.01)

# check covariance is close too
@test 0 < AMP.calcCovarianceBasic(getManifold(X2), getPoints(X2))

@cast pts[i,j] := pts_[j][i]
@test (Statistics.std(pts,dims=2)[1]-1.0) < 3.0
@test (Statistics.std(pts,dims=2)[2]-1.0) < 3.0


##

end

# plotKDE(getBelief(fg, :x2),levels=3)

# spyCliqMat(tree, :x2)

#
