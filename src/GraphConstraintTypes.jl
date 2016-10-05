# this file setup is less than ideal. We will refactor this with packing macros
# to automatically generate the converters and packed data types.
# The eval functions are also currently spread out, they should be concentrated here as well.
# The Pose2D and Pose3D types will most likely be packaged with the RoME package in the future.


# Examples costraint functions which can be used, however,
# \/\/\/these won't pack into ProtoBuf -- need to write special converters
# for Neo4j DataBase storage of complicated types. See constraints hereafter
# for more standard types.
type UniPriorPose2D <: Singleton
  Z::Distributions.MvNormal
end
type GMMPriorPose2D <: Singleton
  Z::Array{Distributions.MvNormal,1}
  W::Array{Float64,1}
end
type KDEPriorPoint2D <: Singleton
  Z::BallTreeDensity
end
type KDERangePoint2D <: Pairwise
  Z::BallTreeDensity
end
# ^^^will work on these later


# Active constraint types listed below
# -------------


# define the pose group
type Odo <: Pairwise
    Zij::Array{Float64,2} # 0rotations, 1translation in each column
    Cov::Array{Float64,2}
    W::Array{Float64,1}
    Odo() = new()
    Odo(x...) = new(x[1], x[2], x[3])
end
type PackedOdo
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
function convert(::Type{FunctionNodeData{PackedOdo}}, d::FunctionNodeData{Odo})
  return FunctionNodeData{PackedOdo}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          convert(PackedOdo, d.fnc))
end
function convert(::Type{FunctionNodeData{Odo}}, d::FunctionNodeData{PackedOdo})
  return FunctionNodeData{Odo}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          convert(Odo, d.fnc))
end
function FNDencode(d::FunctionNodeData{Odo})
  return convert(FunctionNodeData{PackedOdo}, d)
end
function FNDdecode(d::FunctionNodeData{PackedOdo})
  return convert(FunctionNodeData{Odo}, d)
end

type OdoMM <: Pairwise
    Zij::Array{Float64,2} # 0rotations, 1translation in each column
    Cov::Array{Float64,2}
    W::Array{Float64,1}
end


type Ranged <: Pairwise
    Zij::Array{Float64,1}
    Cov::Array{Float64,1}
    W::Array{Float64,1}
    Ranged() = new()
    Ranged(x...) = new(x[1], x[2], x[3])
end
function FNDencode(d::FunctionNodeData{Ranged})
  return d
end
function FNDdecode(d::FunctionNodeData{Ranged})
  return d
end


type GenericMarginal <: Pairwise
    Zij::Array{Float64,1}
    Cov::Array{Float64,1}
    W::Array{Float64,1}
    GenericMarginal() = new()
    GenericMarginal(a,b,c) = new(a,b,c)
end
function FNDencode(d::FunctionNodeData{GenericMarginal})
  return d
end
function FNDdecode(d::FunctionNodeData{GenericMarginal})
  return d
end

# ------------------------------------------------------------


type Obsv2 <: Singleton
    pts::Array{Float64,2}
    bws::Array{Float64,2}
    W::Array{Float64,1}
    Obsv2() = new()
    Obsv2(x...) = new(x[1], x[2], x[3])
end
type PackedObsv2
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
function convert(::Type{FunctionNodeData{PackedObsv2}}, d::FunctionNodeData{Obsv2})
  return FunctionNodeData{PackedObsv2}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          convert(PackedObsv2, d.fnc))
end
function convert(::Type{FunctionNodeData{Obsv2}}, d::FunctionNodeData{PackedObsv2})
  return FunctionNodeData{Obsv2}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          convert(Obsv2, d.fnc))
end
function FNDencode(d::FunctionNodeData{Obsv2})
  return convert(FunctionNodeData{PackedObsv2}, d)
end
function FNDdecode(d::FunctionNodeData{PackedObsv2})
  return convert(FunctionNodeData{Obsv2}, d)
end

# ------------------------------------------------------



type PriorPose2 <: Singleton
    Zi::Array{Float64,2}
    Cov::Array{Float64,2}
    W::Array{Float64,1}
    PriorPose2() = new()
    PriorPose2(x...) = new(x[1], x[2], x[3])
end
type PackedPriorPose2
    vecZij::Array{Float64,1} # 0rotations, 1translation in each column
    dimz::Int64
    vecCov::Array{Float64,1}
    dimc::Int64
    W::Array{Float64,1}
    PackedPriorPose2() = new()
    PackedPriorPose2(x...) = new(x[1], x[2], x[3], x[4], x[5])
end
function convert(::Type{PriorPose2}, d::PackedPriorPose2)
  Zij = reshapeVec2Mat(d.vecZij,d.dimz)
  Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  return PriorPose2(Zij, Cov, d.W)
end
function convert(::Type{PackedPriorPose2}, d::PriorPose2)
  v1 = d.Zij[:];
  v2 = d.Cov[:];
  return PackedPriorPose2(v1,size(d.Zij,1),
                          v2,size(d.Cov,1),
                          d.W)
end
function convert(::Type{FunctionNodeData{PackedPriorPose2}}, d::FunctionNodeData{PriorPose2})
  return FunctionNodeData{PackedPriorPose2}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          convert(PackedPriorPose2, d.fnc))
end
function convert(::Type{FunctionNodeData{PriorPose2}}, d::FunctionNodeData{PackedPriorPose2})
  return FunctionNodeData{PriorPose2}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          convert(PriorPose2, d.fnc))
end
function FNDencode(d::FunctionNodeData{PriorPose2})
  return convert(FunctionNodeData{PackedPriorPose2}, d)
end
function FNDdecode(d::FunctionNodeData{PackedPriorPose2})
  return convert(FunctionNodeData{PriorPose2}, d)
end


# ------------------------------------


type Pose2Pose2 <: Pairwise
    Zij::Array{Float64,2} # 2translations, 1rotation
    Cov::Array{Float64,2}
    W::Array{Float64,1}
    Pose2Pose2() = new()
    Pose2Pose2(x...) = new(x[1], x[2], x[3])
end
type PackedPose2Pose2
  vecZij::Array{Float64,1} # 2translations, 1rotation
  dimz::Int64
  vecCov::Array{Float64,1}
  dimc::Int64
  W::Array{Float64,1}
  PackedPose2Pose2() = new()
  PackedPose2Pose2(x...) = new(x[1], x[2], x[3], x[4], x[5])
end
function convert(::Type{Pose2Pose2}, d::PackedPose2Pose2)
  return Pose2Pose2(reshapeVec2Mat(d.vecZij,d.dimz),
                    reshapeVec2Mat(d.vecCov, d.dimc), W)
end
function convert(::Type{PackedPose2Pose2}, d::Pose2Pose2)
  v1 = d.Zij[:];
  v2 = d.Cov[:];
  return PackedPose2Pose2(v1,size(d.Zij,1),
                          v2,size(d.Cov,1),
                          d.W)
end
function convert(::Type{FunctionNodeData{PackedPose2Pose2}}, d::FunctionNodeData{Pose2Pose2})
  return FunctionNodeData{PackedPose2Pose2}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          convert(PackedPose2Pose2, d.fnc))
end
function convert(::Type{FunctionNodeData{Pose2Pose2}}, d::FunctionNodeData{PackedPose2Pose2})
  return FunctionNodeData{Pose2Pose2}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          convert(Pose2Pose2, d.fnc))
end
function FNDencode(d::FunctionNodeData{Pose2Pose2})
  return convert(FunctionNodeData{PackedPose2Pose2}, d)
end
function FNDdecode(d::FunctionNodeData{PackedPose2Pose2})
  return convert(FunctionNodeData{Pose2Pose2}, d)
end


# --------------------------------------------



type Pose2DPoint2DBearingRange <: Pairwise
    Zij::Array{Float64,2} # bearing and range hypotheses as columns
    Cov::Array{Float64,2}
    W::Array{Float64,1}
    Pose2DPoint2DBearingRange() = new()
    Pose2DPoint2DBearingRange(x...) = new(x[1],x[2],x[3])
end
type PackedPose2DPoint2DBearingRange
    vecZij::Array{Float64,1} # 0rotations, 1translation in each column
    dimz::Int64
    vecCov::Array{Float64,1}
    dimc::Int64
    W::Array{Float64,1}
    PackedPose2DPoint2DBearingRange() = new()
    PackedPose2DPoint2DBearingRange(x...) = new(x[1], x[2], x[3], x[4], x[5])
end
function convert(::Type{Pose2DPoint2DBearingRange}, d::PackedPose2DPoint2DBearingRange)
  Zij = reshapeVec2Mat(d.vecZij,d.dimz)
  Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  return Pose2DPoint2DBearingRange(Zij, Cov, d.W)
end
function convert(::Type{PackedPose2DPoint2DBearingRange}, d::Pose2DPoint2DBearingRange)
  v1 = d.Zij[:];
  v2 = d.Cov[:];
  return PackedPose2DPoint2DBearingRange(v1,size(d.Zij,1),
                                         v2,size(d.Cov,1),
                                         d.W)
end
function convert(::Type{FunctionNodeData{PackedPose2DPoint2DBearingRange}}, d::FunctionNodeData{Pose2DPoint2DBearingRange})
  return FunctionNodeData{PackedPose2DPoint2DBearingRange}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          convert(PackedPose2DPoint2DBearingRange, d.fnc))
end
function convert(::Type{FunctionNodeData{Pose2DPoint2DBearingRange}}, d::FunctionNodeData{PackedPose2DPoint2DBearingRange})
  return FunctionNodeData{Pose2DPoint2DBearingRange}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          convert(Pose2DPoint2DBearingRange, d.fnc))
end
function FNDencode(d::FunctionNodeData{Pose2DPoint2DBearingRange})
  return convert(FunctionNodeData{PackedPose2DPoint2DBearingRange}, d)
end
function FNDdecode(d::FunctionNodeData{PackedPose2DPoint2DBearingRange})
  return convert(FunctionNodeData{Pose2DPoint2DBearingRange}, d)
end


# ------------------------------------------------------
type Pose2DPoint2DRange <: Pairwise
    Zij::Vector{Float64} # bearing and range hypotheses as columns
    Cov::Float64
    W::Vector{Float64}
    Pose2DPoint2DRange() = new()
    Pose2DPoint2DRange(x...) = new(x[1],x[2],x[3])
end
passTypeThrough(d::FunctionNodeData{Pose2DPoint2DRange}) = d



# ------------------------------------------------------
type Point2DPoint2DRange <: Pairwise
    Zij::Vector{Float64} # bearing and range hypotheses as columns
    Cov::Float64
    W::Vector{Float64}
    Point2DPoint2DRange() = new()
    Point2DPoint2DRange(x...) = new(x[1],x[2],x[3])
end
passTypeThrough(d::FunctionNodeData{Point2DPoint2DRange}) = d


# ---------------------------------------------------------

type PriorPoint2D <: Singleton
    mv::MvNormal
    W::Array{Float64,1}
    PriorPoint2D() = new()
    PriorPoint2D(mu, cov, W) = new(MvNormal(mu, cov), W)
end
type PackedPriorPoint2D
    mu::Array{Float64,1}
    vecCov::Array{Float64,1}
    dimc::Int64
    W::Array{Float64,1}
    PackedPriorPoint2D() = new()
    PackedPriorPoint2D(x...) = new(x[1], x[2], x[3], x[4])
end
function convert(::Type{PriorPoint2D}, d::PackedPriorPoint2D)
  Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  return PriorPoint2D(d.mu, Cov, d.W)
end
function convert(::Type{PackedPriorPoint2D}, d::PriorPoint2D)
  v2 = d.mv.Σ.mat[:];
  return PackedPriorPoint2D(d.mv.μ, v2, size(d.mv.Σ.mat,1), d.W)
end

function convert(::Type{FunctionNodeData{PackedPriorPoint2D}}, d::FunctionNodeData{PriorPoint2D})
  return FunctionNodeData{PackedPriorPoint2D}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          convert(PackedPriorPoint2D, d.fnc))
end
function convert(::Type{FunctionNodeData{PriorPoint2D}}, d::FunctionNodeData{PackedPriorPoint2D})
  return FunctionNodeData{PriorPoint2D}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          convert(PriorPoint2D, d.fnc))
end
function FNDencode(d::FunctionNodeData{PriorPoint2D})
  return convert(FunctionNodeData{PackedPriorPoint2D}, d)
end
function FNDdecode(d::FunctionNodeData{PackedPriorPoint2D})
  return convert(FunctionNodeData{PriorPoint2D}, d)
end




# ------------------------------------------------------


type PriorPose3 <: Singleton
    Zi::SE3
    Cov::Array{Float64,2}
    PriorPose3() = new()
    PriorPose3(st::FloatInt, sr::Float64) = new(SE3(0), [[st*eye(3);zeros(3,3)];[zeros(3);sr*eye(3)]])
    PriorPose3(s::SE3, c::Array{Float64,2}) = new(s,c)
end
# type PackedPriorPose3
#     vecZij::Array{Float64,1} # 0rotations, 1translation in each column
#     dimz::Int64
#     vecCov::Array{Float64,1}
#     dimc::Int64
#     W::Array{Float64,1}
#     PackedPriorPose3() = new()
#     PackedPriorPose3(x...) = new(x[1], x[2], x[3], x[4], x[5])
# end
# function convert(::Type{PriorPose3}, d::PackedPriorPose3)
#   Zij = reshapeVec2Mat(d.vecZij,d.dimz)
#   Cov = reshapeVec2Mat(d.vecCov, d.dimc)
#   return PriorPose3(Zij, Cov, d.W)
# end
# function convert(::Type{PackedPriorPose3}, d::PriorPose3)
#   v1 = d.Zij[:];
#   v2 = d.Cov[:];
#   return PackedPriorPose3(v1,size(d.Zij,1),
#                           v2,size(d.Cov,1),
#                           d.W)
# end
# function convert(::Type{FunctionNodeData{PackedPriorPose3}}, d::FunctionNodeData{PriorPose3})
#   return FunctionNodeData{PackedPriorPose3}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           convert(PackedPriorPose3, d.fnc))
# end
# function convert(::Type{FunctionNodeData{PriorPose3}}, d::FunctionNodeData{PackedPriorPose3})
#   return FunctionNodeData{PriorPose3}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           convert(PriorPose3, d.fnc))
# end
# function FNDencode(d::FunctionNodeData{PriorPose3})
#   return convert(FunctionNodeData{PackedPriorPose3}, d)
# end
# function FNDdecode(d::FunctionNodeData{PackedPriorPose3})
#   return convert(FunctionNodeData{PriorPose3}, d)
# end



# ------------------------------------


type Pose3Pose3 <: Pairwise
    Zij::SE3 # 3translations, 3exponential param rotation
    Cov::Array{Float64,2}
    Pose3Pose3() = new()
    Pose3Pose3(st::FloatInt, sr::Float64) = new(SE3(0), [[st*eye(3);zeros(3,3)];[zeros(3);sr*eye(3)]])
    Pose3Pose3(s::SE3, c::Array{Float64,2}) = new(s,c)
end
# type PackedPose3Pose3
#   vecZij::Array{Float64,1} # 2translations, 1rotation
#   dimz::Int64
#   vecCov::Array{Float64,1}
#   dimc::Int64
#   W::Array{Float64,1}
#   PackedPose2Pose2() = new()
#   PackedPose2Pose2(x...) = new(x[1], x[2], x[3], x[4], x[5])
# end
# function convert(::Type{Pose3Pose3}, d::PackedPose3Pose3)
#   return Pose3Pose3(reshapeVec2Mat(d.vecZij,d.dimz),
#                     reshapeVec2Mat(d.vecCov, d.dimc), W)
# end
# function convert(::Type{PackedPose3Pose3}, d::Pose3Pose3)
#   v1 = d.Zij[:];
#   v2 = d.Cov[:];
#   return PackedPose3Pose3(v1,size(d.Zij,1),
#                           v2,size(d.Cov,1),
#                           d.W)
# end
# function convert(::Type{FunctionNodeData{PackedPose3Pose3}}, d::FunctionNodeData{Pose3Pose3})
#   return FunctionNodeData{PackedPose3Pose3}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           convert(PackedPose3Pose3, d.fnc))
# end
# function convert(::Type{FunctionNodeData{Pose3Pose3}}, d::FunctionNodeData{PackedPose3Pose3})
#   return FunctionNodeData{Pose3Pose3}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           convert(Pose3Pose3, d.fnc))
# end
# function FNDencode(d::FunctionNodeData{Pose3Pose3})
#   return convert(FunctionNodeData{PackedPose3Pose3}, d)
# end
# function FNDdecode(d::FunctionNodeData{PackedPose3Pose3})
#   return convert(FunctionNodeData{Pose3Pose3}, d)
# end



# -----------------------
