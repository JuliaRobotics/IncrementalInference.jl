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
  dummy::Int
end

getSample(::Test2, N::Int=1) = (reshape(randn(N),1,N), )

function (tt::Test2)( res::AbstractVector{<:Real},
                      userdata::FactorMetadata,
                      idx::Int,
                      meas::Tuple,
                      tp1::AbstractArray{<:Real,2},
                      tp2::AbstractArray{<:Real,2}  )
  #
  tp1[1,1]=-2.0;
  res[:] .= 1.0
  nothing;
end

##

@testset "CommonConvWrapper test" begin

##

N = 3
tst2 = Test2(1)
A = rand(2,N)
B = rand(2,N)
At = deepcopy(A)
t = Array{Array{Float64,2},1}()
push!(t,A)
push!(t,B)
t[1][1,1] = -10.0
@test A[1,1] == -10

##

fg = initfg()
X0 = addVariable!(fg, :x0, ContinuousEuclid{1})
X1 = addVariable!(fg, :x1, ContinuousEuclid{1})

addFactor!(fg, [:x0;:x1], tst2, graphinit=false)

fmd = FactorMetadata([X0;X1], [:x0; :x1],t, :null, nothing)

##

ccw = CommonConvWrapper(tst2, t[1], 2, t, fmd) # generalwrapper.measurement = rand(1,1)

##

ccw.cpt[Threads.threadid()].factormetadata
ccw.cpt[Threads.threadid()].particleidx
ccw.measurement = (randn(1,N),)
ccw.params[ccw.cpt[Threads.threadid()].activehypo]

x, res = zeros(2), zeros(2)
@time ccw(res, x)


At[1,1] = -2.0
At[2,1] = 0.0
@test A == At

##


end

##

# abstract Nonparametric <: Function
# This is what the internmediate user would be contributing
mutable struct Pose1Pose1Test{T} <: AbstractRelativeRoots
  Dx::T
  Pose1Pose1Test{T}() where T = new()
  Pose1Pose1Test{T}(a::T) where T = new(a)
  # Pose1Pose1Test(a::T) where T = new(a)
end
Pose1Pose1Test(a::T) where T = Pose1Pose1Test{T}(a)

getSample(pp1t::Pose1Pose1Test{T}, N::Int=1) where T = (reshape(rand(pp1t.Dx,N),1,N),)


#proposed standardized parameter list, does not have to be functor
function (Dp::Pose1Pose1Test)(res::AbstractArray{<:Real},
                              userdata::FactorMetadata,
                              idx::Int,
                              meas::Tuple{<:AbstractArray{<:Real,2}},
                              p1::AbstractArray{<:Real,2},
                              p2::AbstractArray{<:Real,2} )
  #
  res[1] = meas[1][1,idx] - (p2[1,idx] - p1[1,idx])
  nothing
end

##

@testset "Test in factor graph setting..." begin

##

N = 100
p1 = rand(1,N)
p2 = rand(1,N)
t = Array{Array{Float64,2},1}()
push!(t,p1)
push!(t,p2)

odo = Pose1Pose1Test(Normal(100.0,1.0))

fg = initfg()
X0 = addVariable!(fg, :x0, ContinuousEuclid{1})
X1 = addVariable!(fg, :x1, ContinuousEuclid{1})
addFactor!(fg, [:x0;:x1], odo, graphinit=false)
fmd = FactorMetadata([X0;X1], [:x0; :x1], t, :null, nothing)

ccw = CommonConvWrapper(odo, t[1], 1, t, fmd, measurement=getSample(odo, N))


freshSamples!(ccw, N, fmd)
x, res = zeros(1), zeros(1)

@time for n in 1:N
  ccw.cpt[Threads.threadid()].particleidx = n

  ccw(res, x)
  # each point should be near 100.0
  @test 50.0 < res[1] < 150.0
end

# which variable position to solve for
ccw.varidx = 2

# do for a single point (i=1)

r = nlsolve( ccw, [0.0], inplace=true )


# for all points
@time for i in 1:N
  ccw.cpt[Threads.threadid()].particleidx = i
    # ccw(x, res)
  r = nlsolve( ccw, ccw.params[ccw.varidx][:,ccw.cpt[Threads.threadid()].particleidx], inplace=true )
  ccw.params[ccw.varidx][1,ccw.cpt[Threads.threadid()].particleidx] = r.zero[1]
end

@test abs(Statistics.mean(p1)-0.0) < 4.0
@test abs(Statistics.mean(p2)-100.0) < 4.0

##

end



@testset "Test with CommonConvWrapper for un-permuted root finding..." begin

##

N = 110
p1 = rand(1,N)
p2 = rand(1,N)
t = Array{Array{Float64,2},1}()
push!(t,p1)
push!(t,p2)

odo = Pose1Pose1Test(Normal(100.0,1.0))
# varidx=2 means we are solving for p2 relative to p1

measurement = getSample(odo, N)
@show zDim = size(measurement[1],1)

solvefor = 2

fg = initfg()
X0 = addVariable!(fg, :x0, ContinuousEuclid{1})
X1 = addVariable!(fg, :x1, ContinuousEuclid{1})
addFactor!(fg, [:x0;:x1], odo, graphinit=false)
fmd = FactorMetadata([X0;X1], [:x0; :x1], t, :null, nothing)

ccw = CommonConvWrapper(odo, t[solvefor], zDim, t, fmd, measurement=measurement)
@show ccw.varidx = solvefor
# gwp = GenericWrapParam{Pose1Pose1Test}(odo, t, 2, 1, (zeros(0,1),) , getSample) #getSample(odo, N)

freshSamples!(ccw, N, fmd)
# and return complete fr/gwp
@time for n in 1:N
  # gwp(x, res)
  ccw.cpt[Threads.threadid()].particleidx = n
  numericSolutionCCW!( ccw )
end

# @show gwp.params

@test 90.0 < Statistics.mean(ccw.params[ccw.varidx]) < 110.0
@test -10.0 < Statistics.mean(ccw.params[1]) < 10.0

println("and in the reverse direction, achieved by simply changing CommonConvWrapper.varidx to 1...")

solvefor = 1
@show ccw.varidx = solvefor
ccw.params[solvefor][:,:] = -100.0*ones(size(ccw.params[solvefor]))
ccw.cpt[Threads.threadid()].X = ccw.params[solvefor]

# @show gwp.params

# fr = FastRootGenericWrapParam{Pose1Pose1Test}(gwp.params[gwp.varidx], zDim, gwp)

freshSamples!(ccw, N, fmd)
@time for n in 1:N
  # gwp(x, res)
  ccw.cpt[Threads.threadid()].particleidx = n
  numericSolutionCCW!( ccw )
end

# @show gwp.params

@test -10.0 < Statistics.mean(ccw.params[1]) < 10.0
@test 90.0 < Statistics.mean(ccw.params[2]) < 110.0

##

end




# use the range only example, should give a circle with nothing in the middle



@testset "Generic convolution testing in factor graph context..." begin

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

tree = resetBuildTree!(fg)

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


end
