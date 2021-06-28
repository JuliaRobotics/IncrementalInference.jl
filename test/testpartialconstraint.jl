# partial constraint development
using Test
using IncrementalInference
# going to introduce two new constraint mutable structs
import IncrementalInference: getSample
using Statistics
using TensorCast


##

mutable struct DevelopPartial{P <: Tuple} <: AbstractPrior
  x::Distribution
  partial::P 
end
getSample(cf::CalcFactor{<:DevelopPartial}, N::Int=1) = (reshape(rand(cf.factor.x, N),1,N), )


mutable struct DevelopDim2 <: AbstractPrior
  x::Distribution
end
getSample(cf::CalcFactor{<:DevelopDim2}, N::Int=1) = (rand(cf.factor.x, N), )


##


N=100 # 50
fg = initfg()


v1 = addVariable!(fg,:x1,ContinuousEuclid{2}(),N=N)

pr = DevelopDim2(MvNormal([0.0;0.0], diagm([0.01;0.01]))) # *Matrix{Float64}(LinearAlgebra.I,2,2)))
f1  = addFactor!(fg,[:x1],pr)

dp = DevelopPartial(Normal(2.0, 1.0),(1,))
f2  = addFactor!(fg,[:x1],dp)


@testset "test evaluation of full constraint prior" begin


pts_ = evalFactor(fg, f1, v1.label, N=N)
@cast pts[i,j] := pts_[j][i]
@test size(pts,1) == 2
@test size(pts,2) == N
@test norm(Statistics.mean(pts,dims=2)[1] .- [0.0]) < 0.3


end

##

memcheck_ = getVal(v1)
@cast memcheck[i,j] := memcheck_[j][i]

@testset "test evaluation of partial constraint prior" begin

X1pts_ = getVal(v1)
@cast X1pts[i,j] := X1pts_[j][i]
pts_ = approxConv(fg, f2.label, :x1, N=N)
@cast pts[i,j] := pts_[j][i]

@test size(pts, 1) == 2
@test size(pts, 2) == N
@test norm(Statistics.mean(pts,dims=2)[1] .- [2.0]) < 0.75
# ensure the correct response from
@test norm(X1pts[1,:] - pts[1,:]) > 2.0
@test norm(X1pts[2,:] - pts[2,:]) < 1e-10
@test norm(X1pts - memcheck) < 1e-10

end

##

@testset "test solving of factor graph" begin


getSolverParams(fg).N = N
tree, smt, hist = solveTree!(fg)
pts_ = getVal(fg, :x1)
@cast pts[i,j] := pts_[j][i]

@test norm(Statistics.mean(pts,dims=2)[1] .- [0.0]) < 0.25
@test norm(Statistics.mean(pts,dims=2)[2] .- [0.0]) < 0.25


end
# plotKDE(getBelief(fg, :x1),levels=3)


##

# @enter predictbelief(fg, :x1, :)


##

mutable struct DevelopPartialPairwise <: AbstractRelativeMinimize
  x::Distribution
  partial::Tuple
  DevelopPartialPairwise(x::Distribution) = new(x, (2,))
end
getSample(cf::CalcFactor{<:DevelopPartialPairwise}, N::Int=1) = (reshape(rand(cf.factor.x, N),1,N), )

function (dp::CalcFactor{<:DevelopPartialPairwise})(meas,
                                                    x1,
                                                    x2  )
  #
  # v0.21+
  return meas[1] - (x2[2]-x1[2])
end



v2 = addVariable!(fg,:x2,ContinuousEuclid{2},N=N)


dpp = DevelopPartialPairwise(Normal(10.0, 1.0))
f3  = addFactor!(fg,[:x1;:x2],dpp)


dp2 = DevelopPartial( Normal(-20.0, 1.0), (1,) )
f4  = addFactor!(fg,[:x2;], dp2, graphinit=false)
doautoinit!(fg, :x2)

# drawGraph(fg, show=true)


##

@testset "test evaluation of multiple simultaneous partial constraints" begin
global fg

##

initAll!(fg)
valx2_ = getVal(fg, :x2)
@cast valx2[i,j] := valx2_[j][i]
pts_ = approxConv(fg, :x1x2f1, :x2, N=N) # evalFactor(fg, f3, v2.index, N=N)
@cast pts[i,j] := pts_[j][i]
@test size(pts,1) == 2
@test norm(Statistics.mean(pts,dims=2)[2] .- [10.0]) < 3.0
@test norm(valx2[1,:] - pts[1,:]) < 1e-5

pts_ = approxConv(fg, :x2f1, :x2, N=N) # evalFactor(fg, f4, v2.index, N=N)
@cast pts[i,j] := pts_[j][i]
@test size(pts,1) == 2
@test norm(Statistics.mean(pts,dims=2)[1] .- [-20.0]) < 0.75
@test (Statistics.std(pts,dims=2)[1] .- 1.0) < 0.4

##

end

##

# keep previous values to ensure funciton evaluation is modifying correct data fields

@warn "restore findRelatedFromPotential as testset!"
# @testset "test findRelatedFromPotential..." begin
# global v2, fg, f3, f4, N


thefac = getFactor(fg, :x1x2f1)

X2lpts_ = getVal(getVariable(fg, :x2))
@cast X2lpts[i,j] := X2lpts_[j][i]
keepaside, = findRelatedFromPotential(fg, thefac, :x2, N=N)
@test Ndim(keepaside) == 2
lpts_ = getPoints(keepaside)
@cast lpts[i,j] := lpts_[j][i]
@test length(lpts_) == N

@show X2lpts[2,95:100]
@show lpts[2,95:100]
@show getPoints(keepaside)

# DevelopPartialPairwise must only modify the second dimension of proposal distribution on X2
@test norm(X2lpts[1,:] - lpts[1,:]) < 1e-10
# @test norm(X2lpts[2,:] - lpts[2,:]) > 1e-10 # 10*N*0.5 # why so big?
memcheck_ = getVal(v2)
@cast memcheck[i,j] := memcheck_[j][i]
@test norm(X2lpts - memcheck) < 1e-10


X2lpts_ = getVal(v2)
@cast X2lpts[i,j] := X2lpts_[j][i]
p4, = findRelatedFromPotential(fg, f4, v2.label, N=N)
@test Ndim(p4) == 2
lpts_ = getPoints(keepaside)
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


# combination of partials
val_, = predictbelief(fg, v2, [f3;f4], N=N)
@cast val[i,j] := val_[j][i]
# plotKDE(kde!(val),levels=3)
@test norm(Statistics.mean(val,dims=2)[1] .- [-20.0]) < 3.0
@test norm(Statistics.mean(val,dims=2)[2] .- [10.0]) < 3.0
@test (Statistics.std(val,dims=2)[1] .- 1.0) < 3.0
@test (Statistics.std(val,dims=2)[2] .- 1.0) < 3.0



getSolverParams(fg).N = N
tree, smt, hist = solveTree!(fg)


pts_ = getVal(fg, :x1)
@cast pts[i,j] := pts_[j][i]
@test norm(Statistics.mean(pts,dims=2)[1] .- [0.0]) < 0.5
@test norm(Statistics.mean(pts,dims=2)[2] .- [0.0]) < 0.5

pts_ = getVal(fg, :x2)
@cast pts[i,j] := pts_[j][i]
@test norm(Statistics.mean(pts,dims=2)[1] .- [-20.0]) < 3.0
@test norm(Statistics.mean(pts,dims=2)[2] .- [10.0]) < 3.0
@test (Statistics.std(pts,dims=2)[1]-1.0) < 3.0
@test (Statistics.std(pts,dims=2)[2]-1.0) < 3.0

##

end

# plotKDE(getBelief(fg, :x2),levels=3)

# spyCliqMat(tree, :x2)

#
