# this file setup is less than ideal. We will refactor this with packing macros
# to automatically generate the converters and packed data structs.
# The eval functions are also currently spread out, they should be concentrated here as well.
# The Pose2D and Pose3D structs will most likely be packaged with the RoME package in the future.




# Active constraint structs listed below
# -------------


# define the simple 1D odo
# TODO -- rework to use Distributions rather than Z and Cov
"""
$(TYPEDEF)
"""
mutable struct Odo <: FunctorPairwise
    Zij::Array{Float64,2} # 0rotations, 1translation in each column
    Cov::Array{Float64,2}
    W::Array{Float64,1}
    Odo() = new()
    Odo(x...) = new(x[1], x[2], x[3])
end
# TODO -- only computing first node
function (odo::Odo)(res::Vector{Float64},
    userdata::FactorMetadata,
    idx::Int,
    meas::Tuple{Array{Float64,2}},
    p1::Array{Float64},
    p2::Array{Float64}  )

  res[1] = meas[1][1,idx] - (p2[1,idx] - p1[1,idx])
  nothing
end
function getSample(odo::Odo, N::Int=1)
  (reshape(rand(Distributions.Normal(odo.Zij[1,1], odo.Cov[1,1]), N ),1,N),)
end
# function getSample(odo::Odo, N::Int=1)
#   ret = zeros(1,N)
#   if size(odo.Zij,2) > 1
#     error("getSample(::Odo,::Int) can't handle multi-column at present")
#   end
#   for i in 1:N
#     ret[1,i] = odo.Cov[1]*randn()+odo.Zij[1]
#   end
#   # rand(Distributions.Normal(odo.Zij[1],odo.Cov[1]), N)'
#   return ret
# end
"""
$(TYPEDEF)
"""
mutable struct PackedOdo <: PackedInferenceType
    vecZij::Array{Float64,1} # 0rotations, 1translation in each column
    dimz::Int
    vecCov::Array{Float64,1}
    dimc::Int
    W::Array{Float64,1}
    PackedOdo() = new()
    PackedOdo(x...) = new(x[1], x[2], x[3], x[4], x[5])
end
# function convert(::Type{Odo}, d::PackedOdo)
#   Zij = reshapeVec2Mat(d.vecZij,d.dimz)
#   Cov = reshapeVec2Mat(d.vecCov, d.dimc)
#   return Odo(Zij, Cov, d.W)
# end
# function convert(::Type{PackedOdo}, d::Odo)
#   v1 = d.Zij[:];
#   v2 = d.Cov[:];
#   return PackedOdo(v1,length(v1),
#                     v2,length(v2),
#                     d.W)
# end



"""
$(TYPEDEF)
"""
mutable struct OdoMM <: Pairwise
    Zij::Array{Float64,2} # 0rotations, 1translation in each column
    Cov::Array{Float64,2}
    W::Array{Float64,1}
end
function getSample(odo::OdoMM, N::Int=1)
  ret = zeros(1,N)
  if size(odo.Zij,2) > 1
    error("getSample(::OdoMM,::Int) can't handle multi-column at present")
  end
  for i in 1:N
    ret[1,i] = odo.Cov[1]*randn()+odo.Zij[1]
  end
  # rand(Distributions.Normal(odo.Zij[1],odo.Cov[1]), N)'
  return (ret,)
end


"""
$(TYPEDEF)
"""
mutable struct Ranged <: FunctorPairwise
    Zij::Array{Float64,1}
    Cov::Array{Float64,1}
    W::Array{Float64,1}
    Ranged() = new()
    Ranged(x...) = new(x[1], x[2], x[3])
end
"""
$(TYPEDEF)
"""
mutable struct PackedRanged <: PackedInferenceType
    Zij::Array{Float64,1}
    Cov::Array{Float64,1}
    W::Array{Float64,1}
    PackedRanged() = new()
    PackedRanged(x...) = new(x[1], x[2], x[3])
end
function convert(::Type{Ranged}, r::PackedRanged)
  return Ranged(r.Zij, r.Cov, r.W)
end
function convert(::Type{PackedRanged}, r::Ranged)
  return PackedRanged(r.Zij, r.Cov, r.W)
end
function (ra::Ranged)(res::Vector{Float64},
    userdata::FactorMetadata,
    idx::Int,
    meas::Tuple{Array{Float64,2}},
    p1::Array{Float64},
    l1::Array{Float64})

  res[1] = meas[1][1,idx] - abs(l1[1,idx] - p1[1,idx])
  nothing
end
function getSample(ra::Ranged, N::Int=1)
  ret = zeros(1,N)
  for i in 1:N
    ret[1,i] = ra.Cov[1]*randn()+ra.Zij[1]
  end
  # rand(Distributions.Normal(odo.Zij[1],odo.Cov[1]), N)'
  return (ret,)
end


"""
$(TYPEDEF)
"""
mutable struct GenericMarginal <: FunctorPairwise
    Zij::Array{Float64,1}
    Cov::Array{Float64,1}
    W::Array{Float64,1}
    GenericMarginal() = new()
    GenericMarginal(a,b,c) = new(a,b,c)
end
"""
$(TYPEDEF)
"""
mutable struct PackedGenericMarginal <: PackedInferenceType
    Zij::Array{Float64,1}
    Cov::Array{Float64,1}
    W::Array{Float64,1}
    PackedGenericMarginal() = new()
    PackedGenericMarginal(a,b,c) = new(a,b,c)
end
function convert(::Type{PackedGenericMarginal}, d::GenericMarginal)
  return PackedGenericMarginal(d.Zij, d.Cov, d.W)
end
function convert(::Type{GenericMarginal}, d::PackedGenericMarginal)
  return GenericMarginal(d.Zij, d.Cov, d.W)
end

# ------------------------------------------------------------
