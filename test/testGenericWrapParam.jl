# test GenericWrapParam
using Base: Test
using Distributions
using NLsolve
using KernelDensityEstimate
using IncrementalInference

import IncrementalInference: getSample

println("FunctorWorks")
type FunctorWorks
  a::Array{Float64,2}
end

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



println("FunctorArray")

type FunctorArray{T}
  # This was a tuple, but array will like work better in the long term
  fnc!::Function
  a::Array{T, 1}
end

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




println("GenericWrapParam test")


function test2(res::Vector{Float64}, idx::Int, meas::Tuple{Array{Float64,2}}, tp1::Array{Float64,2}, tp2::Array{Float64,2})
  tp1[1,1]=-2.0;
  res[:] = 1.0
  nothing;
end


A = rand(2,3)
B = rand(2,3)
At = deepcopy(A)
t = Array{Array{Float64,2},1}()
push!(t,A)
push!(t,B)
# @show typeof(t)
generalwrapper = GenericWrapParam{typeof(test2)}(test2, t, 1, 1)
# generalwrapper.measurement = rand(1,1)

x, res = zeros(2), zeros(2)
@time generalwrapper(x,res)

At[1,1] = -2.0
At[2,1] = 0.0
@test A == At





println("Test in factor graph setting...")

# abstract Nonparametric <: Function
# This is what the internmediate user would be contributing
type Pose1Pose1Test{T} <: FunctorPairwise
  Dx::T
  Pose1Pose1Test() = new()
  # Pose1Pose1Test{T}(::Int) = new(T())
  Pose1Pose1Test{T}(a::T) = new(a)
end
getSample{T}(pp1t::Pose1Pose1Test{T}, N::Int=1) = (rand(pp1t.Dx,N)',)

#proposed standardized parameter list, does not have to be functor
function (Dp::Pose1Pose1Test)(res::Array{Float64},
      idx::Int,
      meas::Tuple{Array{Float64,2}},
      p1::Array{Float64,2},
      p2::Array{Float64,2} )
  #
  # println("Dp::Pose1Pose1Test, in-place idx=$(idx)")
  res[1] = meas[1][1,idx] - (p2[1,idx] - p1[1,idx])
  nothing
end

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
  generalwrapper(x, res)
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


# function evalPotential(factor::GenericWrapParam, Xi::Array{Graphs.ExVertex,1}, solveforid::Int64; N:Int=100)
#
#
# end

println("Test with FastRootGenericWrapParam for un-permuted root finding...")
function _tmp()
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

# ARR = Array{Array{Float64,2},1}()
# maxlen, sfidx = prepareparamsarray!(ARR, Xi, N, solvefor)
@show gwp.varidx
gwp.measurement = gwp.samplerfnc(gwp.usrfnc!, N)
@show zDim = size(gwp.measurement,1)
fr = FastRootGenericWrapParam{Pose1Pose1Test}(gwp.params[gwp.varidx], zDim, gwp)
# and return complete fr/gwp

@time for gwp.particleidx in 1:N
  # gwp(x, res)
  numericRootGenericRandomizedFnc!( fr )
end

function _fooSeq(fr)
  for i in 1:N
    fr.gwp.particleidx = i
  # gwp(x, res)
    numericRootGenericRandomizedFnc!( fr )
  end
end

function _fooPar(fr)
  Threads.@threads for i in 1:N
  ccall(:jl_, Void, (Any,), i)
  fr.gwp.particleidx = i
  # gwp(x, res)
    numericRootGenericRandomizedFnc!( fr )
  end
end

_fooSeq(fr)
@time _fooSeq(fr)
_fooPar(fr)
@time _fooPar(fr)

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


println("Test with FastRootGenericWrapParam for permuted root finding...")
warn("not implemented yet")

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

v1=addNode!(fg, :x1, p1, N=N)
v2=addNode!(fg, :x2, p2, N=N)
f1 = addFactor!(fg, [v1], Obsv2(p1, getBW(d1)[:,1]', [1.0]))

odo = Pose1Pose1Test{Normal}(Normal(100.0,1.0))
f2 = addFactor!(fg, [v1;v2], odo)

tree = wipeBuildNewTree!(fg)

pts = getVal(getVert(fg,:x1))
@test abs(Base.mean(pts)-0.0) < 10.0
pts = getVal(getVert(fg,:x2))
@test abs(Base.mean(pts)-0.0) < 10.0

@time [inferOverTreeR!(fg, tree) for i in 1:3]


# using Gadfly
# plot(y=rand(10))
# plotKDE(getVertKDE(fg,:x2))

pts = getVal(getVert(fg,:x1))
@test abs(Base.mean(pts)-0.0) < 10.0

pts = getVal(getVert(fg,:x2))
@test abs(Base.mean(pts)-100.0) < 10.0





# repeat tests with SolverUtilities version
