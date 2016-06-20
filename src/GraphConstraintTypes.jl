# this file setup is less than ideal. We will refactor this with packing macros
# to automatically generate the converters and packed data types.
# The eval functions are also currently spread out, they should be concentrated here as well.
# The Pose2D and Pose3D types will most likely be packaged with the RoME package in the future.


# \/\/\/these won't pack into ProtoBuf, serialize might work
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


# -------------


# define the pose group
type Odo <: Pairwise
    Zij::Array{Float64,2} # 0rotations, 1translation in each column
    Cov::Array{Float64,2}
    W::Array{Float64,1}
end
type PackedOdo
    vecZij::Array{Float64,1} # 0rotations, 1translation in each column
    dimz::Int64
    vecCov::Array{Float64,1}
    dimc::Int64
    W::Array{Float64,1}
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
  return FunctionNodeData{PackedOdo}(d.fncargvID, d.eliminated, d.potentialused,
          convert(PackedOdo, d.fnc))
end
function convert(::Type{FunctionNodeData{Odo}}, d::FunctionNodeData{PackedOdo})
  return FunctionNodeData{Odo}(d.fncargvID, d.eliminated, d.potentialused,
          convert(Odo, d.fnc))
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
end

type GenericMarginal <: Pairwise
    Zij::Array{Float64,1}
    Cov::Array{Float64,1}
    W::Array{Float64,1}
    GenericMarginal() = new()
    GenericMarginal(a,b,c,d) = new(a,b,c,d)
end


# ------------------------------------------------------------


type Obsv2 <: Singleton
    pts::Array{Float64,2}
    bws::Array{Float64,2}
    W::Array{Float64,1}
end
type PackedObsv2
    vecZij::Array{Float64,1} # 0rotations, 1translation in each column
    dimz::Int64
    vecCov::Array{Float64,1}
    dimc::Int64
    W::Array{Float64,1}
end
function convert(::Type{Obsv2}, d::PackedObsv2)
  Zij = reshapeVec2Mat(d.vecZij,d.dimz)
  Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  return Obsv2(Zij, Cov, d.W)
end
function convert(::Type{PackedObsv2}, d::Obsv2)
  v1 = d.pts[:];
  v2 = d.bws[:];
  return PackedObsv2(v1,length(v1),
                    v2,length(v2),
                    d.W)
end
function convert(::Type{FunctionNodeData{PackedObsv2}}, d::FunctionNodeData{Obsv2})
  return FunctionNodeData{PackedObsv2}(d.fncargvID, d.eliminated, d.potentialused,
          convert(PackedObsv2, d.fnc))
end
function convert(::Type{FunctionNodeData{Obsv2}}, d::FunctionNodeData{PackedObsv2})
  return FunctionNodeData{Obsv2}(d.fncargvID, d.eliminated, d.potentialused,
          convert(Obsv2, d.fnc))
end


# ------------------------------------------------------



type PriorPose2 <: Singleton
    Zi::Array{Float64,2}
    Cov::Array{Float64,2}
    W::Array{Float64,1}
end
type PackedPriorPose2
    vecZij::Array{Float64,1} # 0rotations, 1translation in each column
    dimz::Int64
    vecCov::Array{Float64,1}
    dimc::Int64
    W::Array{Float64,1}
end
function convert(::Type{PriorPose2}, d::PackedPriorPose2)
  Zij = reshapeVec2Mat(d.vecZij,d.dimz)
  Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  return PriorPose2(Zij, Cov, d.W)
end
function convert(::Type{PackedPriorPose2}, d::PriorPose2)
  v1 = d.Zij[:];
  v2 = d.Cov[:];
  return PackedPriorPose2(v1,length(v1),
                    v2,length(v2),
                    d.W)
end
function convert(::Type{FunctionNodeData{PackedPriorPose2}}, d::FunctionNodeData{PriorPose2})
  return FunctionNodeData{PackedPriorPose2}(d.fncargvID, d.eliminated, d.potentialused,
          convert(PackedPriorPose2, d.fnc))
end
function convert(::Type{FunctionNodeData{PriorPose2}}, d::FunctionNodeData{PackedPriorPose2})
  return FunctionNodeData{PriorPose2}(d.fncargvID, d.eliminated, d.potentialused,
          convert(PriorPose2, d.fnc))
end


# ------------------------------------




type Pose2Pose2 <: Pairwise
    Zij::Array{Float64,2} # 2translations, 1rotation
    Cov::Array{Float64,2}
    W::Array{Float64,1}
end
type PackedPose2Pose2
  vecZij::Array{Float64,1} # 2translations, 1rotation
  dimz::Int64
  vecCov::Array{Float64,1}
  dimc::Int64
  W::Array{Float64,1}
end
function convert(::Type{Pose2Pose2}, d::PackedPose2Pose2)
  return Pose2Pose2(reshapeVec2Mat(d.vecZij,d.dimz),
                    reshapeVec2Mat(d.vecCov, d.dimc), W)
end
function convert(::Type{PackedPose2Pose2}, d::Pose2Pose2)
  v1 = d.vecZij[:];
  v2 = d.Cov[:];
  return PackedPose2Pose2(v1,length(v1),
                          v2,length(v2),
                          d.W)
end
function convert(::Type{FunctionNodeData{PackedPose2Pose2}}, d::FunctionNodeData{Pose2Pose2})
  return FunctionNodeData{PackedPose2Pose2}(d.fncargvID, d.eliminated, d.potentialused,
          convert(PackedPose2Pose2, d.fnc))
end
function convert(::Type{FunctionNodeData{Pose2Pose2}}, d::FunctionNodeData{PackedPose2Pose2})
  return FunctionNodeData{Pose2Pose2}(d.fncargvID, d.eliminated, d.potentialused,
          convert(Pose2Pose2, d.fnc))
end



# --------------------------------------------



type Pose2DPoint2DBearingRange <: Pairwise
    Zij::Array{Float64,2} # bearing and range hypotheses as columns
    Cov::Array{Float64,2}
    W::Array{Float64,1}
end
type PackedPose2DPoint2DBearingRange
    vecZij::Array{Float64,1} # 0rotations, 1translation in each column
    dimz::Int64
    vecCov::Array{Float64,1}
    dimc::Int64
    W::Array{Float64,1}
end
function convert(::Type{Pose2DPoint2DBearingRange}, d::PackedPose2DPoint2DBearingRange)
  Zij = reshapeVec2Mat(d.vecZij,d.dimz)
  Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  return Pose2DPoint2DBearingRange(Zij, Cov, d.W)
end
function convert(::Type{PackedPose2DPoint2DBearingRange}, d::Pose2DPoint2DBearingRange)
  v1 = d.Zij[:];
  v2 = d.Cov[:];
  return PackedPose2DPoint2DBearingRange(v1,length(v1),
                                        v2,length(v2),
                                        d.W)
end
function convert(::Type{FunctionNodeData{PackedPose2DPoint2DBearingRange}}, d::FunctionNodeData{Pose2DPoint2DBearingRange})
  return FunctionNodeData{PackedPose2DPoint2DBearingRange}(d.fncargvID, d.eliminated, d.potentialused,
          convert(PackedPose2DPoint2DBearingRange, d.fnc))
end
function convert(::Type{FunctionNodeData{Pose2DPoint2DBearingRange}}, d::FunctionNodeData{PackedPose2DPoint2DBearingRange})
  return FunctionNodeData{Pose2DPoint2DBearingRange}(d.fncargvID, d.eliminated, d.potentialused,
          convert(Pose2DPoint2DBearingRange, d.fnc))
end



# -----------------------
