# test numericRootGenericRandomized
using IncrementalInference
using Test

import IncrementalInference: getSample

##

mutable struct LineResidual <: AbstractRelativeRoots
  m::Float64
  c::Float64
end

getSample(lr::LineResidual, N::Int=1) = (reshape(rand(Normal(),N),1,N),)

# y = mx + c
# res = y      -      m*x  -  c
#       meas          variables and fixed values
function (lr::LineResidual)(res::AbstractVector{<:Real},
                            fmd::FactorMetadata,
                            idx::Int,
                            z::Tuple,
                            x::AbstractArray{<:Real,2},
                            y::AbstractArray{<:Real,2}  )
  #
  res[1] = z[1][idx] - (y[1,idx] - (lr.m*x[1,idx] + lr.c))
  nothing
end

##


@testset "test CommonConvWrapper{T}, solve of residual functions..." begin

##

function assembleConvType(dfgvrs::Vector{<:DFGVariable}, functor::T, xDim::Int, zDim::Int, nvars::Int) where {T <: AbstractRelativeRoots}
  # @info "assembleConvType -- development testing function only, not intended for production."
  N = 3

  vars = Array{Array{Float64,2},1}()
  for i in 1:nvars
    push!(vars, zeros(xDim, N))
  end

  fmd = IIF._defaultFactorMetadata(dfgvrs, arrRef=vars)

  CommonConvWrapper(functor,vars[1],zDim,vars, fmd, measurement=(zeros(zDim,N),))
end


##

# TODO -- expand testing to include subcomponent tests from numericSolutionCCW!
lr1 = LineResidual(2.0, 3.0)

fg = initfg()
X0 = addVariable!(fg, :x0, ContinuousScalar)
X1 = addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0,:x1], lr1, graphinit=false)

ccw = assembleConvType([X0;X1], lr1, 1, 1, 2)

res = zeros(1)
ccw(res, zeros(1))

# gwp(x, res)
ccw.cpt[1].particleidx = 1
numericSolutionCCW!( ccw )

@test abs(ccw.cpt[1].Y[1] + 1.50) < 1e-5


##

end


# println("Test shuffling function")
#
# function testshuffle!(res, x)
#   # println("testshuffle!(x,res) gets x=$(x), and length(res)=$(length(res))")
#   @show res
#   @show x
#   nval = sum(abs.(x-[11.0;12.0;13.0]))
#   # trick the solver so that f/f' = 0/1e10, and therefore < tol
#   # resulting in termination of solver
#   if nval < 1e-8
#     res[1:2] = [0.0;0.0] # fake at root evaluation
#   else
#     res[1:2] = 1e10*randn(2) # fake large gradient
#   end
#   nothing
# end
#
# # test shuffling function
# ALLTESTS = Bool[]
# @show x0 = collect(11:13)+0.0
# for i in 1:10
#   y = numericRootGenericRandomizedFnc(
#           testshuffle!,
#           2, 3, x0,
#           perturb=0.0    )
#   @show y
#   @show y.%10
#   @show yy = y.%10.0
#   push!(ALLTESTS, norm(yy-collect(1:3)) < 0.1)
# end
# @test sum(ALLTESTS) > 9
#
#
# println("Test if shuffling results in correct mapping for solving")
#
#
# function testshuffle2!(res, x)
#   # println("testshuffle2!(x,res) gets x=$(x), and length(res)=$(length(res))")
#   res[1:2] = x - [1.0;2.0] # fake at root evaluation
#   nothing
# end
#
# # test shuffling function
# x0 = collect(11:12)+0.0
# for i in 1:10
#   println("starting with x0=$(x0)")
#   y = numericRootGenericRandomizedFnc(
#           testshuffle2!,
#           2, 2, x0,
#           perturb=0.0,
#           testshuffle=true    )
#   println("result y=$(y)")
#   @test norm(y-collect(1:2)) < 0.1
# end
#
#
#
# # x is dimension 3, z dimension is 2
# function residualrandtest!(res, x)
#   val = norm(x[1:2])
#   res[1] = 10.0 - val
#   nothing
# end
#
# for i in 1:10
#   x0 = ones(2)
#   y = numericRootGenericRandomizedFnc(
#           residualrandtest!,
#           1, 2, x0   )
#   #
#   @test abs(norm(y) - 10.0) < 1e-4
# end
