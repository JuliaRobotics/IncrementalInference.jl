# this file setup is less than ideal. We will refactor this with packing macros
# to automatically generate the converters and packed data types.
# The eval functions are also currently spread out, they should be concentrated here as well.
# The Pose2D and Pose3D types will most likely be packaged with the RoME package in the future.




# Active constraint types listed below
# -------------


# define the pose group
type Odo <: FunctorPairwise
    Zij::Array{Float64,2} # 0rotations, 1translation in each column
    Cov::Array{Float64,2}
    W::Array{Float64,1}
    Odo() = new()
    Odo(x...) = new(x[1], x[2], x[3])
end
function (odo::Odo)(res::Vector{Float64},
    idx::Int,
    meas::Tuple{Array{Float64,2}},
    p1::Array{Float64},
    p2::Array{Float64}  )

  res[1] = meas[1][1,idx] - (p2[1,idx] - p1[1,idx])
  nothing
end
function getSample(odo::Odo, N::Int=1)
  (rand(Distributions.Normal(0.0, odo.Cov[1,1]), N )',)
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
type PackedOdo <: PackedInferenceType
    vecZij::Array{Float64,1} # 0rotations, 1translation in each column
    dimz::Int64
    vecCov::Array{Float64,1}
    dimc::Int64
    W::Array{Float64,1}
    PackedOdo() = new()
    PackedOdo(x...) = new(x[1], x[2], x[3], x[4], x[5])
end
function convert(::Type{Odo}, d::PackedOdo)
  Zij = reshapeVec2Mat(d.vecZij,d.dimz)
  Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  return Odo(Zij, Cov, d.W)
end
function convert(::Type{PackedOdo}, d::Odo)
  v1 = d.Zij[:];
  v2 = d.Cov[:];
  return PackedOdo(v1,length(v1),
                    v2,length(v2),
                    d.W)
end



type OdoMM <: Pairwise
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


type Ranged <: FunctorPairwise
    Zij::Array{Float64,1}
    Cov::Array{Float64,1}
    W::Array{Float64,1}
    Ranged() = new()
    Ranged(x...) = new(x[1], x[2], x[3])
end
type PackedRanged <: PackedInferenceType
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


type GenericMarginal <: FunctorPairwise
    Zij::Array{Float64,1}
    Cov::Array{Float64,1}
    W::Array{Float64,1}
    GenericMarginal() = new()
    GenericMarginal(a,b,c) = new(a,b,c)
end
type PackedGenericMarginal <: PackedInferenceType
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


type Obsv2 <: FunctorSingleton
    pts::Array{Float64,2}
    bws::Array{Float64,2}
    W::Array{Float64,1}
    Obsv2() = new()
    Obsv2(x...) = new(x[1], x[2], x[3])
end
type PackedObsv2 <: PackedInferenceType
    vecZij::Array{Float64,1} # 0rotations, 1translation in each column
    dimz::Int64
    vecCov::Array{Float64,1}
    dimc::Int64
    W::Array{Float64,1}
    PackedObsv2() = new()
    PackedObsv2(x...) = new(x[1],x[2],x[3],x[4],x[5])
end
function convert(::Type{Obsv2}, d::PackedObsv2)
  Zij = reshapeVec2Mat(d.vecZij,d.dimz)
  Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  return Obsv2(Zij, Cov, d.W)
end
function convert(::Type{PackedObsv2}, d::Obsv2)
  v1 = d.pts[:];
  v2 = d.bws[:];
  return PackedObsv2(v1,size(d.pts,1),
                    v2,size(d.bws,1),
                    d.W)
end
function getSample(z::Obsv2, N::Int=1)
  return (KernelDensityEstimate.sample(kde!(z.pts, z.bws), N),)
end
