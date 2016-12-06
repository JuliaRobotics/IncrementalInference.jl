# test numericRootGenericRandomized
using IncrementalInference, TransformUtils
using Base.Test

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

y = numericRootGenericRandomized(lineResidual!, 1, [1.0], [2.0;3.0], randn(1))
@test abs(y[1] + 1.0) < 1e-10

y = numericRootGenericRandomized(lineResidual!, 1, [3.0], [2.0;3.0], randn(1))
@test abs(y[1] - 0.0) < 1e-10

y = numericRootGenericRandomized(lineResidual!, 1, [5.0], [2.0;3.0], randn(1))
@test abs(y[1] - 1.0) < 1e-10




println("Increased dimension test")

# 3 dimensional line, z = [a b][x y]' + c
function rotationResidual!(res::Vector{Float64}, z::Vector{Float64}, var::Tuple)
  # fixed = var[1]
  # state = var[2]
  q1 = convert(Quaternion, Euler(z...))
  q2 = convert(Quaternion, so3(var[2]))
  qq = q1*q_conj(q2)
  res[1:3] = vee(convert(so3, qq))
  nothing
end

# 3 dimensional line, z - ([ax ay][x y]' + c)
for i in 1:10
  eul = 0.25*randn(3)
  y = numericRootGenericRandomized(rotationResidual!, 3, eul, zeros(3), 0.1*randn(3))
  # test the result
  @test TransformUtils.compare(convert(Quaternion, Euler(eul...)),
                              convert(Quaternion, so3(y)), tol=1e-8)
end

for i in 1:10
  eul = 0.25*randn(3)
  gg = (x, res) -> rotationResidual!(res, eul, (zeros(0),x))
  y = numericRootGenericRandomizedFnc(
          gg,
          3, 3, 0.1*randn(3)    )
  # test the result
  @test TransformUtils.compare(convert(Quaternion, Euler(eul...)),
                              convert(Quaternion, so3(y)), tol=1e-8)
end

# type CameraIntrinsic
#   K::Array{Float64,2}
#   CameraIntrinsic(::Void) = new()
#   CameraIntrinsic(;x0=320.0,y0=240.0,fx=510.0,fy=510.0,s=0.0) = new([[fx;s;x0]';[0.0;fy;y0]';[0.0;0;1]'])
# end
#
# # Camera extrinsic must be world in camera frame (cRw)
# type CameraExtrinsic
#   R::SO3
#   t::Vector{Float64}
#   CameraExtrinsic(::Void) = new()
#   CameraExtrinsic(;R=SO3(0),t=zeros(3)) = new(R, t)
# end
# type CameraModelFull
#   ci::CameraIntrinsic
#   ce::CameraExtrinsic
#   # cd::CameraDistortion
#   CameraModelFull(::Void) = new()
#   CameraModelFull(;ci=CameraIntrinsic(), ce=CameraExtrinsic()) = new(ci,ce)
# end
# function project!(ret::Vector{Float64}, ci::CameraIntrinsic, ce::CameraExtrinsic, pt::Vector{Float64})
#   res = ci.K*(ce.R.R*pt + ce.t)
#   ret[1:2] = res[1:2]./res[3]
#   nothing
# end
# project!(ret::Vector{Float64}, cm::CameraModelFull, pt::Vector{Float64}) = project!(ret, cm.ci, cm.ce, pt)
# function project(cm::CameraModelFull, pt::Vector{Float64})
#   res = Vector{Float64}(2)
#   project!(res, cm, pt)
#   return res
# end
#
# # pinhole camera model
# # (x, y)/f = (X, Y)/Z
# function cameraResidual!(
#       res::Vector{Float64},
#       z::Vector{Float64},
#       ci::CameraIntrinsic,
#       ce::CameraExtrinsic,
#       pt::Vector{Float64}  )
#   # in place memory operations
#   project!(res, ci, ce, pt)
#   res[1:2] .*= -1.0
#   res[1:2] += z[1:2]
#   nothing
# end
#
# ci = CameraIntrinsic()
# ce = CameraExtrinsic()
# pt = [1.0;0.0;5.0]
#
#
# gg = (x, res) -> cameraResidual!(res, x, ci, ce, pt)
# # res = zeros(2)
# # @time gg([0.0;0.0], res)
#
# # Profile.clear()
# # @profile
# y = numericRootGenericRandomizedFnc(
#         gg,
#         2, 2, randn(2)  )
# #
# @test abs((ci.K[1,3]+ci.K[1,1]*pt[1]/pt[3]) - y[1]) < 1e-8
# @test abs(y[2] - ci.K[2,3]) < 1e-8

# using ProfileView
# ProfileView.view()




# test convexity of rotation residual
#
# N=100
# X = linspace(-pi/2,pi/2,N)
# eul = 0.5*randn(3)
# # RES = zeros(3, N)
# # for i in 1:N
# #   res = zeros(3)
# #   rotationResidual!(res, eul, (zeros(0), [X[i];0.0;0.0]) )
# #   RES[:,i] = res
# # end
#
# using Gadfly
#
# function easyplot(x,y, eul)
#   res = zeros(3)
#   rotationResidual!(res, eul, (zeros(0), [x;0.0;y]) )
#   return norm(res)
# end
#
# plot(z = (x,y) -> easyplot(x,y,eul), x=X, y=X, Geom.contour)
