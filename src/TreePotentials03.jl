# 3 dimensionsal evaluation functions.
# Transforms code is scattered, and will be moved to better location with stable 3D Examples
# also worth it to consolidate if other affine transforms package is available

# using Euler angles for linear sampling in product of potentials
# Gaussian model for prior
function evalPotential(obs::PriorPose3, Xi::Array{Graphs.ExVertex,1}; N::Int64=200)
  mu = veeEuler(obs.Zi)
  return rand( MvNormal(mu, obs.Cov), N)
end

# Project all particles (columns) Xval with Z, that is for all  SE3(Xval[:,i])*Z
function projectParticles(Xval::Array{Float64,2}, Z::SE3, Cov::Array{Float64,2})
  # TODO optimize for more speed with threads and better memory management
  r,c = size(Xval)
  RES = zeros(r,c) #*cz

  ent, x = SE3(0), SE3(0)
  ENT = rand( MvNormal(zeros(6), Cov), c )

  j=1
  # for j in 1:cz
    for i in 1:c
      x.R, x.t = convert(SO3,Euler(Xval[4:6,i])), Xval[1:3,i]
      ent.R, ent.t = convert(SO3,so3(ENT[4:6,i])), ENT[1:3,i]
      newval = x*Z*ent
      RES[1:r,i*j] = veeEuler(newval)
    end
  # end
  #
  return RES
end

# Still limited to linear sampler, then reprojected onto ball -- TODO upgrade manifold sampler
function evalPotential(odom::Pose3Pose3, Xi::Array{Graphs.ExVertex,1}, Xid::Int64; N::Int64=200)
    # rz,cz = size(odom.Zij)
    Xval = Array{Float64,2}()
    # implicit equation portion -- bi-directional pairwise function made explicit here
    if Xid == Xi[1].index #odom.
        # reverse direction
        println(" reverse direction")
        Z = inverse(odom.Zij)
        Xval = getVal(Xi[2])
    elseif Xid == Xi[2].index
        # forward direction
        println(" forward direction")
        Z = odom.Zij
        Xval = getVal(Xi[1])
    else
        error("Bad evalPairwise Pose3Pose3")
    end

    return projectParticles(Xval, Z, odom.Cov)
end
