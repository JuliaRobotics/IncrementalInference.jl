# test GenericWrapParam
using Base: Test
using Distributions
using NLsolve
using KernelDensityEstimate
using IncrementalInference

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


struct Test2 <: FunctorInferenceType
end

@testset "GenericWrapParam test" begin


function (tt::Test2)(res::Vector{Float64}, userdata::FactorMetadata, idx::Int, meas::Tuple{Array{Float64,2}}, tp1::Array{Float64,2}, tp2::Array{Float64,2})
  tp1[1,1]=-2.0;
  res[:] = 1.0
  nothing;
end

tst2 = Test2()
A = rand(2,3)
B = rand(2,3)
At = deepcopy(A)
t = Array{Array{Float64,2},1}()
push!(t,A)
push!(t,B)
t[1][1,1] = -10.0
@test A[1,1] == -10
# @show typeof(t)
generalwrapper = GenericWrapParam{Test2}(tst2, t, 1, 1)
ccw = CommonConvWrapper(tst2, t[1], 2, t) # generalwrapper.measurement = rand(1,1)


x, res = zeros(2), zeros(2)
@show t[1][:,1]
@time generalwrapper(res, x)

At[1,1] = -2.0
At[2,1] = 0.0
@test A == At


x, res = zeros(2), zeros(2)
@time ccw(res, x)


At[1,1] = -2.0
At[2,1] = 0.0
@test A == At

# resource requirements

# fm = FactorMetadata()
# meas = (randn(2,3), )
#
# @time tst2(res, fm, 1, meas, A, B)

end


# abstract Nonparametric <: Function
# This is what the internmediate user would be contributing
mutable struct Pose1Pose1Test{T} <: FunctorPairwise
  Dx::T
  Pose1Pose1Test{T}() where T = new()
  Pose1Pose1Test{T}(a::T) where T = new(a)
  # Pose1Pose1Test(a::T) where T = new(a)
end
Pose1Pose1Test(a::T) where T = Pose1Pose1Test{T}(a)

getSample{T}(pp1t::Pose1Pose1Test{T}, N::Int=1) = (reshape(rand(pp1t.Dx,N),1,N),)


#proposed standardized parameter list, does not have to be functor
function (Dp::Pose1Pose1Test)(res::Array{Float64},
            userdata::FactorMetadata,
            idx::Int,
            meas::Tuple{Array{Float64,2}},
            p1::Array{Float64,2},
            p2::Array{Float64,2} )
  #
  res[1] = meas[1][1,idx] - (p2[1,idx] - p1[1,idx])
  nothing
end


@testset "Test in factor graph setting..." begin


N = 101
p1 = rand(1,N)
p2 = rand(1,N)
t = Array{Array{Float64,2},1}()
push!(t,p1)
push!(t,p2)

odo = Pose1Pose1Test{Normal}(Normal(100.0,1.0))
generalwrapper = GenericWrapParam{Pose1Pose1Test}(odo, t, 1, 1, getSample(odo, N), getSample)
# generalwrapper.measurement = getSample(odo, N)

x, res = zeros(1), zeros(1)
@time for generalwrapper.particleidx in 1:N
    # nlsolve( generalwrapper, x0? .. )
  generalwrapper(res, x)
  # each point should be near 100.0
  @test res[1] > 50.0
end


println("Test with NLsolve for root finding using generalwrapper functor.")

generalwrapper.varidx = 2
@time for i in 1:N
  generalwrapper.particleidx = i
    # generalwrapper(x, res)
  r = nlsolve( generalwrapper, generalwrapper.params[generalwrapper.varidx][:,generalwrapper.particleidx] )
  generalwrapper.params[generalwrapper.varidx][1,generalwrapper.particleidx] = r.zero[1]
end

@test abs(Base.mean(p1)-0.0) < 3.0
@test abs(Base.mean(p2)-100.0) < 3.0



## Repeat for new CommonConvWrapper

N = 101
p1 = rand(1,N)
p2 = rand(1,N)
t = Array{Array{Float64,2},1}()
push!(t,p1)
push!(t,p2)

odo = Pose1Pose1Test{Normal}(Normal(100.0,1.0))
ccw = CommonConvWrapper(odo, t[1], 1, t, measurement=getSample(odo, N))

x, res = zeros(1), zeros(1)

@time for ccw.particleidx in 1:N
    # nlsolve( generalwrapper, x0? .. )
  ccw(res, x)
  # each point should be near 100.0
  @test res[1] > 50.0
end

ccw.varidx = 2
@time for i in 4:N
  ccw.particleidx = i
    # ccw(x, res)
  r = nlsolve( ccw, ccw.params[ccw.varidx][:,ccw.particleidx] )
  ccw.params[ccw.varidx][1,ccw.particleidx] = r.zero[1]
end

@test abs(Base.mean(p1)-0.0) < 3.0
@test abs(Base.mean(p2)-100.0) < 3.0

end


# function evalPotential(factor::GenericWrapParam, Xi::Array{Graphs.ExVertex,1}, solveforid::Int; N:Int=100)
#
#
# end

@testset "Test with FastRootGenericWrapParam for un-permuted root finding..." begin

N = 110
p1 = rand(1,N)
p2 = rand(1,N)
t = Array{Array{Float64,2},1}()
push!(t,p1)
push!(t,p2)

odo = Pose1Pose1Test{Normal}(Normal(100.0,1.0))
# varidx=2 means we are solving for p2 relative to p1
gwp = GenericWrapParam{Pose1Pose1Test}(odo, t, 2, 1, (zeros(0,1),) , getSample) #getSample(odo, N)

# fgr = FastGenericRoot{GenericWrapParam{Pose1Pose1Test}}(1, 1, generalwrapper)
# numericRootGenericRandomizedFnc!( fgr )
# # but this only represents the last step, should happen internal -- TODO
#
# @show p2

@show gwp.varidx
gwp.measurement = gwp.samplerfnc(gwp.usrfnc!, N)
@show zDim = size(gwp.measurement,1)
fr = FastRootGenericWrapParam{Pose1Pose1Test}(gwp.params[gwp.varidx], zDim, gwp)

# and return complete fr/gwp
@time for gwp.particleidx in 1:N
  # gwp(x, res)
  numericRootGenericRandomizedFnc!( fr )
end

# @show gwp.params

@test 90.0 < Base.mean(gwp.params[gwp.varidx]) < 110.0
@test -10.0 < Base.mean(gwp.params[1]) < 10.0

println("and in the reverse direction, achieved by simply changing GenericWrapParam.varidx to 1...")

@show gwp.varidx = 1
gwp.params[1][:,:] = -100.0*ones(size(gwp.params[1]))

# @show gwp.params

fr = FastRootGenericWrapParam{Pose1Pose1Test}(gwp.params[gwp.varidx], zDim, gwp)

@time for gwp.particleidx in 1:100
  # gwp(x, res)
  numericRootGenericRandomizedFnc!( fr )
end

# @show gwp.params

@test -10.0 < Base.mean(gwp.params[1]) < 10.0
@test 90.0 < Base.mean(gwp.params[2]) < 110.0

end


@testset "Test with FastRootGenericWrapParam for permuted root finding..." begin

warn("test not implemented yet")

end

# use the range only example, should give a circle with nothing in the middle



println("GenericWrapParam testing in factor graph context...")

N=101
p1 = randn(1,N)
d1 = kde!(p1)
p2 = randn(1,N)
t = Array{Array{Float64,2},1}()
push!(t,p1)
push!(t,p2)

fg = emptyFactorGraph()
# fg.registeredModuleFunctions[:Main] = getSample

v1=addNode!(fg, :x1, ContinuousScalar, N=N)
v2=addNode!(fg, :x2, ContinuousScalar, N=N)
bws = getBW(d1)[:,1]
f1 = addFactor!(fg, [v1], Obsv2(p1, reshape(bws, 1, length(bws)), [1.0]))

odo = Pose1Pose1Test{Normal}(Normal(100.0,1.0))
f2 = addFactor!(fg, [v1;v2], odo)

tree = wipeBuildNewTree!(fg)

pts = getVal(getVert(fg,:x1))
@test abs(Base.mean(pts)-0.0) < 10.0
pts = getVal(getVert(fg,:x2))
@test abs(Base.mean(pts)-0.0) < 10.0

inferOverTreeR!(fg, tree, N=N)
inferOverTreeR!(fg, tree)
# @time [inferOverTreeR!(fg, tree, N=N) for i in 1:3];


# using Gadfly
# plot(y=rand(10))
# plotKDE(getVertKDE(fg,:x2))

pts = getVal(getVert(fg,:x1))
@test abs(Base.mean(pts)-0.0) < 10.0

pts = getVal(getVert(fg,:x2))
@test abs(Base.mean(pts)-100.0) < 10.0





# repeat tests with SolverUtilities version
