# this file setup is less than ideal. We will refactor this with packing macros
# to automatically generate the converters and packed data types.
# The eval functions are also currently spread out, they should be concentrated here as well.
# The Pose2D and Pose3D types will most likely be packaged with the RoME package in the future.


# heavy use of multiple dispatch for converting between packed and original data types during DB usage
function convert{T <: InferenceType, P <: PackedInferenceType}(::Type{FunctionNodeData{T}}, d::PackedFunctionNodeData{P})
  return FunctionNodeData{T}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          Symbol(d.frommodule), convert(T, d.fnc))
end
function convert{P <: PackedInferenceType, T <: InferenceType}(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T})
  return PackedFunctionNodeData{P}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          string(d.frommodule), convert(P, d.fnc))
end

function FNDencode{T <: InferenceType, P <: PackedInferenceType}(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T})
  return convert(PackedFunctionNodeData{P}, d) #PackedFunctionNodeData{P}
end
function FNDdecode{T <: InferenceType, P <: PackedInferenceType}(::Type{FunctionNodeData{T}}, d::PackedFunctionNodeData{P})
  return convert(FunctionNodeData{T}, d) #FunctionNodeData{T}
end

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


type Ranged <: Pairwise
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



type GenericMarginal <: Pairwise
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


type Obsv2 <: Singleton
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
