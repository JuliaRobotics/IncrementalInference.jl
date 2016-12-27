# test GenericWrapParam
using Base: Test
using IncrementalInference
using Distributions


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



println("FunctorTuple")

# but this does not work
type FunctorTuple
  fnc!::Function
  a::Tuple
end

function testtuple!(a1::Array{Float64, 2}, a2::Array{Float64,2})
  a1[1,1] = -1.0
  @show a1
  nothing
end

function (fi::FunctorTuple)(x)
  fi.fnc!(fi.a...)
end

A = rand(2,3)
B = rand(2,3)
t = (A,B)
At = deepcopy(A)
At[1,1] = -1.0

fvar = FunctorTuple(testtuple!, t)
fvar([0.0])
@test At == A


# this would be in SolverUtilities
# this will likely expand with more internal bells and whistles
# to perform in place memory operations for array values in
# type GenericWrapParam <: Function
#   usrfnc!::Function
#   params::Tuple
#   varidx::Int
#   particleidx::Int
# end

# potential functor approach
# function (p::GenericWrapParam)(x, res)
#   # approximates by not considering cross indices among parameters
#   p.params[p.varidx][:, p.particleidx] = x
#   p.usrfnc!(res, p.particleidx, p.params...)
# end



println("GenericWrapParam test")
# quick test
A = randn(2, 3)
B = randn(2, 3)

At = deepcopy(A)

# this would get called

function test2(res::Vector{Float64}, idx::Int, tp1::Array{Float64,2}, tp2::Array{Float64,2})
  tp1[1,1]=-2.0;
  res[:] = 1.0
  nothing;
end

t = (A,B)
# @show typeof(t)
generalwrapper = GenericWrapParam(test2, t, 1, 1)

v = zeros(2)
@time generalwrapper(v,v)

At[1,1] = -2.0
At[2,1] = 0.0
#
@test A == At



# abstract Nonparametric <: Function
# This is what the internmediate user would be contributing
type Pose1Pose1Test{T} <: Function
  Dx::T
  Pose1Pose1Test() = new()
  Pose1Pose1Test{T}(::Int) = new(T())
  Pose1Pose1Test{T}(a::T) = new(a)
end

#proposed standardized form functor interaction
function (Dp::Pose1Pose1Test)(res::Array{Float64},
      idx::Int,
      p1::Array{Float64,2},
      p2::Array{Float64,2} )
  #
  res[1] = rand(Dp.Dx) - (p2[1,idx] - p1[1,idx])

  nothing
end

odo = Pose1Pose1Test{Normal}(Normal())
generalwrapper = GenericWrapParam(test2, t, 1, 1)



# repeat tests with SolverUtilities version
