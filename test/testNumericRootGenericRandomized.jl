# test numericRootGenericRandomized
using IncrementalInference
using Base.Test

# TODO -- resid function definitions are changing
# residual functions must have the form:
# function residual!(
      # residual::Vector{Float64},
      # z::Vector{Float64},
      # variables::Tuple )



# y = mx + c
# res = y      -      m*x  -  c
#       meas          variables and fixed values
function lineResidual!(res::Vector{Float64}, y::Vector{Float64}, var::Tuple)
  m, c = var[1][1], var[1][2]
  x = var[2]
  res[1] = y[1] - (m*x[1] + c)
  nothing
end

# and calling
# function numericRootGenericRandomized(
      # residFnc::Function,
      # zDim::Int,
      # measurement::Vector{Float64},
      # fixed::Vector{Float64},
      # x0::Vector{Float64};
      # perturb::Float64=0.01  )

params = [2.0;3.0]
gg = (x, res) -> lineResidual!(res, [1.0], (params, x))
y = numericRootGenericRandomizedFnc(gg, 1, 1, randn(1))
# y = numericRootGenericRandomized(lineResidual!, 1, [1.0], params, randn(1))
@test abs(y[1] + 1.0) < 1e-10

gg = (x, res) -> lineResidual!(res, [3.0], (params, x))
y = numericRootGenericRandomizedFnc(gg, 1, 1, randn(1))
# y = numericRootGenericRandomized(lineResidual!, 1, [3.0], params, randn(1))
@test abs(y[1] - 0.0) < 1e-10

gg = (x, res) -> lineResidual!(res, [5.0], (params, x))
y = numericRootGenericRandomizedFnc(gg, 1, 1, randn(1))
# y = numericRootGenericRandomized(lineResidual!, 1, [5.0], params, randn(1))
@test abs(y[1] - 1.0) < 1e-10




# println("Increased dimension test")
#
# # 3 dimensional line, z = [a b][x y]' + c
# function rotationresidual!(res::Vector{Float64}, z::Vector{Float64}, var::Tuple)
#   q1 = convert(Quaternion, Euler(z...))
#   q2 = convert(Quaternion, so3(var[2]))
#   qq = q1*q_conj(q2)
#   res[1:3] = vee(convert(so3, qq))
#   nothing
# end
#
# for i in 1:10
#   eul = 0.25*randn(3)
#   gg = (x, res) -> rotationresidual!(res, eul, (zeros(0),x))
#   y = numericRootGenericRandomizedFnc(
#           gg,
#           3, 3, 0.1*randn(3)    )
#   # test the result
#   @test TransformUtils.compare(convert(Quaternion, Euler(eul...)),
#                               convert(Quaternion, so3(y)), tol=1e-8)
# end

println("Test shuffling function")

function testshuffle!(x, res)
  # println("testshuffle!(x,res) gets x=$(x), and length(res)=$(length(res))")
  nval = sum(abs(x-[11.0;12.0;13.0]))
  # trick the solver so that f/f' = 0/1e10, and therefore < tol
  # resulting in termination of solver
  if nval < 1e-8
    res[1:2] = [0.0;0.0] # fake at root evaluation
  else
    res[1:2] = 1e10*randn(2) # fake large gradient
  end
  nothing
end

# test shuffling function
ALLTESTS = Bool[]
x0 = collect(11:13)+0.0
for i in 1:10
  y = numericRootGenericRandomizedFnc(
          testshuffle!,
          2, 3, x0,
          perturb=0.0    )
  @show y
  @show y.%10
  @show yy = y.%10.0
  push!(ALLTESTS, norm(yy-collect(1:3)) < 0.1)
end
@test sum(ALLTESTS) > 9


println("Test if shuffling results in correct mapping for solving")


function testshuffle2!(x, res)
  # println("testshuffle2!(x,res) gets x=$(x), and length(res)=$(length(res))")
  res[1:2] = x - [1.0;2.0] # fake at root evaluation
  nothing
end

# test shuffling function
x0 = collect(11:12)+0.0
for i in 1:10
  println("starting with x0=$(x0)")
  y = numericRootGenericRandomizedFnc(
          testshuffle2!,
          2, 2, x0,
          perturb=0.0,
          testshuffle=true    )
  println("result y=$(y)")
  @test norm(y-collect(1:2)) < 0.1
end



# x is dimension 3, z dimension is 2
function residualrandtest!(x, res)
  val = norm(x[1:2])
  res[1] = 10.0 - val
  nothing
end

for i in 1:10
  x0 = ones(2)
  y = numericRootGenericRandomizedFnc(
          residualrandtest!,
          1, 2, x0   )
  #
  @test abs(norm(y) - 10.0) < 1e-4
end



warn("Test FastRootGenericWrapParam{T} not implemented yet")





# test convexity of rotation residual
#
# N=100
# X = linspace(-pi/2,pi/2,N)
# eul = 0.5*randn(3)
# # RES = zeros(3, N)
# # for i in 1:N
# #   res = zeros(3)
# #   rotationresidual!(res, eul, (zeros(0), [X[i];0.0;0.0]) )
# #   RES[:,i] = res
# # end
#
# using Gadfly
#
# function easyplot(x,y, eul)
#   res = zeros(3)
#   rotationresidual!(res, eul, (zeros(0), [x;0.0;y]) )
#   return norm(res)
# end
#
# plot(z = (x,y) -> easyplot(x,y,eul), x=X, y=X, Geom.contour)
