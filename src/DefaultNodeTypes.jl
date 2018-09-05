# default node types in IIF

SamplableBelief = Union{Distributions.Distribution, KernelDensityEstimate.BallTreeDensity, AliasingScalarSampler}


struct ContinuousScalar <: InferenceVariable
  dims::Int
  labels::Vector{String}
  ContinuousScalar() = new(1, String["";])
end
struct ContinuousMultivariate <: InferenceVariable
  dims::Int
  labels::Vector{String}
  ContinuousMultivariate() = new()
  ContinuousMultivariate(x) = new(x, String["";])
end


struct Prior{T} <: IncrementalInference.FunctorSingleton where T <: SamplableBelief
  Z::T
end
getSample(s::Prior, N::Int=1) = (reshape(rand(s.Z,N),:,N), )

struct LinearConditional{T} <: IncrementalInference.FunctorPairwise where T <: SamplableBelief
  Z::T
end
getSample(s::LinearConditional, N::Int=1) = (reshape(rand(s.Z,N),:,N), )
function (s::LinearConditional)(res::Array{Float64},
                                userdata::FactorMetadata,
                                idx::Int,
                                meas::Tuple{Array{Float64, 2}},
                                X1::Array{Float64,2},
                                X2::Array{Float64,2}  )
  #
  res[1] = meas[1][idx] - (X2[1,idx] - X1[1,idx])
  nothing
end

struct MixtureLinearConditional <: IncrementalInference.FunctorPairwise
  Z::Vector{Distribution}
  C::Categorical
end
getSample(s::MixtureLinearConditional, N::Int=1) = (reshape.(rand.(s.Z, N),1,:)..., rand(s.C, N))
function (s::MixtureLinearConditional)(res::Array{Float64},
                               userdata::FactorMetadata,
                               idx::Int,
                               meas::Tuple,
                               X1::Array{Float64,2},
                               X2::Array{Float64,2}  )
  #
  res[1] = meas[meas[end][idx]][idx] - (X2[1,idx] - X1[1,idx])
  nothing
end



## packed types are still developed by hand.  Future versions would likely use a @packable macro to write Protobuf safe versions of factors

struct PackedPrior <: PackedInferenceType
  Z::String
  PackedPrior() = new()
  PackedPrior(z::AS) where {AS <: AbstractString} = new(z)
end
function convert(::Type{PackedPrior}, d::Prior)
  PackedPrior(string(d.Z))
end
function convert(::Type{Prior}, d::PackedPrior)
  Prior(extractdistribution(d.Z))
end


struct PackedLinearConditional <: PackedInferenceType
  Z::String
  PackedLinearConditional() = new()
  PackedLinearConditional(z::AS) where {AS <: AbstractString} = new(z)
end
function convert(::Type{PackedLinearConditional}, d::LinearConditional)
  PackedLinearConditional(string(d.Z))
end
function convert(::Type{LinearConditional}, d::PackedLinearConditional)
  LinearConditional(extractdistribution(d.Z))
end


struct PackedMixtureLinearConditional <: PackedInferenceType
  strs::Vector{String}
  cat::String
  PackedMixtureLinearConditional() = new()
  PackedMixtureLinearConditional(z::Vector{<:AbstractString}, cstr::AS) where {AS <: AbstractString} = new(z, cstr)
end
function convert(::Type{PackedMixtureLinearConditional}, d::MixtureLinearConditional)
  PackedMixtureLinearConditional(string.(d.Z), string(d.C))
end
function convert(::Type{MixtureLinearConditional}, d::PackedMixtureLinearConditional)
  MixtureLinearConditional(extractdistribution.(d.strs), extractdistribution(d.cat))
end
