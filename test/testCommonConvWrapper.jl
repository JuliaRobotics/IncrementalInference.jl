# test CommonConvWrapper

using Test
# using NLsolve
using IncrementalInference
using Manifolds
using Statistics
using TensorCast

import IncrementalInference: getSample, getManifold

mutable struct FunctorWorks
  a::Array{Float64,2}
end

@testset "FunctorWorks" begin

function (fw::FunctorWorks)(x)
  fw.a[1,1] = -1.0
  nothing
end

A = rand(2,3)
At = deepcopy(A)
At[1,1] = -1.0

fvar = FunctorWorks(A)
fvar(0.0)
@test At == A

end

mutable struct FunctorArray{T}
  # This was a tuple, but array will like work better in the long term
  fnc!::Function
  a::Array{T, 1}
end

##

@testset "FunctorArray" begin

##

function testarray!(a1::Array{Float64, 2}, a2::Array{Float64,2})
  a1[1,1] = -1.0
  @show a1
  nothing
end

function (fi::FunctorArray)(x)
  fi.fnc!(fi.a...)
end

A = rand(2,3)
B = rand(2,3)
t = Array{Array{Float64,2},1}()
push!(t,A)
push!(t,B)
At = deepcopy(A)
At[1,1] = -1.0

fvar = FunctorArray(testarray!, t)
fvar([0.0])
@test At == A

##

end

##

# abstract Nonparametric <: Function
# This is what the internmediate user would be contributing
mutable struct Pose1Pose1Test{T} <: AbstractManifoldMinimize
  Dx::T
end

getManifold(::Pose1Pose1Test) = TranslationGroup(1)

function getSample(cf::CalcFactor{<:Pose1Pose1Test})
  return rand(cf.factor.Dx, 1)
end

#proposed standardized parameter list, does not have to be functor
function (cf::CalcFactor{<:Pose1Pose1Test})(Dx,
                                            p1,
                                            p2 )
  #
  return Dx[1] - (p2[1] - p1[1])
end

##


@testset "Test with CommonConvWrapper for un-permuted root finding..." begin

##

N = 110
p1 = rand(1,N)
p2 = rand(1,N)
t = Array{Array{Float64,2},1}()
push!(t,p1)
push!(t,p2)

odo = Pose1Pose1Test(Normal(100.0,1.0))


fg = initfg()
X0 = addVariable!(fg, :x0, ContinuousEuclid{1})
initVariable!(fg, :x0, [zeros(1) for _ in 1:100])
X1 = addVariable!(fg, :x1, ContinuousEuclid{1})
addFactor!(fg, [:x0;:x1], odo, graphinit=false)


pts = approxConv(fg, getFactor(fg, :x0x1f1), :x1)

## Now check the contents of internal CCW

ccw = IIF._getCCW(fg, :x0x1f1)

ptr_ = ccw.varValsAll[][ccw.varidx[]]
@cast tp1[i,j] := ptr_[j][i]
@test 90.0 < Statistics.mean(tp1) < 110.0
ptr_ = ccw.varValsAll[][1]
@cast tp2[i,j] := ptr_[j][i]
@test -10.0 < Statistics.mean(tp2) < 10.0

##

println("and in the reverse direction")

initVariable!(fg, :x1, [100*ones(1) for _ in 1:100])

pts = approxConv(fg, getFactor(fg, :x0x1f1), :x0)

ptr_ = ccw.varValsAll[][1]
@cast tp1[i,j] := ptr_[j][i]
@test -10.0 < Statistics.mean(tp1) < 10.0
ptr_ = ccw.varValsAll[][2]
@cast tp2[i,j] := ptr_[j][i]
@test 90.0 < Statistics.mean(tp2) < 110.0

##

end


# use the range only example, should give a circle with nothing in the middle


@testset "Generic convolution testing in factor graph context..." begin

##

N=100
p1 = [randn(1) for _ in 1:N]
d1 = manikde!(TranslationGroup(1), p1)
p2 = [randn(1) for _ in 1:N]
t = Vector{Vector{Vector{Float64}}}()
push!(t,p1)
push!(t,p2)

fg = initfg()

v1=addVariable!(fg, :x1, ContinuousScalar, N=N)
v2=addVariable!(fg, :x2, ContinuousScalar, N=N)
bws = getBW(d1)[:,1]
f1 = addFactor!(fg, [v1], Prior(manikde!(TranslationGroup(1), p1, bw=bws)) )

odo = Pose1Pose1Test(Normal(100.0,1.0))
f2 = addFactor!(fg, [v1;v2], odo)

tree = buildTreeReset!(fg)

pts_ = getBelief(fg,:x1) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test abs(Statistics.mean(pts)-0.0) < 10.0
pts_ = getBelief(fg,:x2) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test abs(Statistics.mean(pts)-0.0) < 10.0

##

tree = solveTree!(fg)

##

pts_ = getBelief(fg,:x1) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test abs(Statistics.mean(pts)-0.0) < 10.0

pts_ = getBelief(fg,:x2) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test abs(Statistics.mean(pts)-100.0) < 10.0

##

end
