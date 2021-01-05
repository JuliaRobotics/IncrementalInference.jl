# test CommonConvWrapper

using Test
using NLsolve
using IncrementalInference
using Statistics

import IncrementalInference: getSample

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

end


##

# abstract Nonparametric <: Function
# This is what the internmediate user would be contributing
mutable struct Pose1Pose1Test{T} <: AbstractRelativeRoots
  Dx::T
  # Pose1Pose1Test{T}() where T = new()
  # Pose1Pose1Test{T}(a::T) where T = new(a)
  # Pose1Pose1Test(a::T) where T = new(a)
end
# Pose1Pose1Test(a::T) where T = Pose1Pose1Test{T}(a)

getSample(cf::CalcFactor{<:Pose1Pose1Test}, N::Int=1) = (reshape(rand(cf.factor.Dx,N),1,N),)


#proposed standardized parameter list, does not have to be functor
function (cf::CalcFactor{<:Pose1Pose1Test})(res::AbstractVector{<:Real},
                                            Dx,
                                            p1,
                                            p2 )
  #
  res[1] = Dx[1] - (p2[1] - p1[1])
  nothing
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
initManual!(fg, :x0, zeros(1,100))
X1 = addVariable!(fg, :x1, ContinuousEuclid{1})
addFactor!(fg, [:x0;:x1], odo, graphinit=false)


pts = approxConv(fg, getFactor(fg, :x0x1f1), :x1)

## Now check the contents of internal CCW

ccw = IIF._getCCW(fg, :x0x1f1)

@test 90.0 < Statistics.mean(ccw.params[ccw.varidx]) < 110.0
@test -10.0 < Statistics.mean(ccw.params[1]) < 10.0

##

println("and in the reverse direction")

initManual!(fg, :x1, 100*ones(1,100))

pts = approxConv(fg, getFactor(fg, :x0x1f1), :x0)

@test -10.0 < Statistics.mean(ccw.params[1]) < 10.0
@test 90.0 < Statistics.mean(ccw.params[2]) < 110.0

##

end


# use the range only example, should give a circle with nothing in the middle


@testset "Generic convolution testing in factor graph context..." begin

##

N=100
p1 = randn(1,N)
d1 = kde!(p1)
p2 = randn(1,N)
t = Array{Array{Float64,2},1}()
push!(t,p1)
push!(t,p2)

fg = initfg()
# fg.registeredModuleFunctions[:Main] = getSample

v1=addVariable!(fg, :x1, ContinuousScalar, N=N)
v2=addVariable!(fg, :x2, ContinuousScalar, N=N)
bws = getBW(d1)[:,1]
f1 = addFactor!(fg, [v1], Prior(kde!(p1, bws)) )

odo = Pose1Pose1Test(Normal(100.0,1.0))
f2 = addFactor!(fg, [v1;v2], odo)

tree = buildTreeReset!(fg)

pts = getVal(getVariable(fg,:x1))
@test abs(Statistics.mean(pts)-0.0) < 10.0
pts = getVal(getVariable(fg,:x2))
@test abs(Statistics.mean(pts)-0.0) < 10.0

# inferOverTreeR!(fg, tree, N=N)
# inferOverTreeR!(fg, tree)
# @time [inferOverTreeR!(fg, tree, N=N) for i in 1:3];
tree, smt, hist = solveTree!(fg)


# using Gadfly
# plot(y=rand(10))
# plotKDE(getBelief(fg,:x2))

pts = getVal(getVariable(fg,:x1))
@test abs(Statistics.mean(pts)-0.0) < 10.0

pts = getVal(getVariable(fg,:x2))
@test abs(Statistics.mean(pts)-100.0) < 10.0

##

end
