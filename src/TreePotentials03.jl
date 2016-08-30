# 3 dimensionsal evaluation functions.
# Transforms code is scattered, and will be moved to better location with stable 3D Examples
# also worth it to consolidate if other affine transforms package is available

function unpackSE3AxisAngle(s::SE3)
  aa = convert(AxisAngle, s.R)
  v = zeros(7)
  v[1:3] = s.t
  # theta, axis
  v[4] = aa.theta
  v[5:7] = aa.ax[:]
  return v
end


# DX = [transx, transy, theta]
function addPose3Pose3(x::Array{Float64,1}, dx::Array{Float64,1})
    X = SE3(x)
    DX = SE3(dx)
    return unpackSE3AxisAngle(X*DX)
end


# Still limited to linear sampler, then reprojected onto ball -- TODO upgrade manifold sampler
function evalPotential(odom::Pose3Pose3, Xi::Array{Graphs.ExVertex,1}, Xid::Int64)
    error("Not implemented yet")
    # rz,cz = size(odom.Zij)
    # Xval = Array{Float64,2}()
    # # implicit equation portion -- bi-directional pairwise function
    # if Xid == Xi[1].index #odom.
    #     Z = se3vee(SE3(vec(odom.Zij)) \ eye(4))
    #     Xval = getVal(Xi[2])
    # elseif Xid == Xi[2].index
    #     Z = odom.Zij
    #     Xval = getVal(Xi[1])
    # else
    #     error("Bad evalPairwise Pose2Pose2")
    # end
    #
    # r,c = size(Xval)
    # RES = zeros(r,c*cz)
    #
    # # TODO -- this should be the covariate error from Distributions, only using diagonals here (approxmition for speed in first implementation)
    # ENT = randn(r,c)
    # for d in 1:r
    #    @fastmath @inbounds ENT[d,:] = ENT[d,:].*odom.Cov[d,d]
    # end
    # # Repeat ENT values for new modes from meas
    # for j in 1:cz
    #   for i in 1:c
    #       z = Z[1:r,j].+ENT[1:r,i]
    #       RES[1:r,i*j] = addPose3Pose3(Xval[1:r,i], z )
    #   end
    # end
    #
    # return RES
end
