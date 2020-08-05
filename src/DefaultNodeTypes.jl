# default node types in IIF

SamplableBelief = Union{Distributions.Distribution, KernelDensityEstimate.BallTreeDensity, AliasingScalarSampler}

#Supported types for parametric
ParametricTypes = Union{Normal, MvNormal}


"""
$(TYPEDEF)

Most basic continuous scalar variable in a `::DFG.AbstractDFG` object.
"""
struct ContinuousScalar <: InferenceVariable
  dims::Int
  manifolds::Tuple{Symbol}
  # ContinuousScalar(;manifolds::Tuple{Symbol}=(:Euclid,)) = new(1, manifolds)
  function ContinuousScalar(;manifolds=nothing)
    manifolds !== nothing &&
      Base.depwarn("ContinuousScalar keyword argument manifolds is deprecated.", :ContinuousScalar)

    return new(1, (:Euclid,))
  end
end

"""
$(TYPEDEF)

Continuous variable of dimension `.dims` on manifold `.manifolds`.
"""
struct ContinuousMultivariate{T1 <: Tuple} <: InferenceVariable
  dims::Int
  manifolds::T1
  ContinuousMultivariate{T}() where {T} = new()
  ContinuousMultivariate{T}(x::Int; manifolds::T=(:Euclid,)) where {T <: Tuple} = new(x, manifolds)
end

function ContinuousMultivariate(x::Int;
                                manifolds::T1=(:Euclid,)  )  where {T1 <: Tuple}
  #
  maniT = length(manifolds) < x ? ([manifolds[1] for i in 1:x]...,) : manifolds
  ContinuousMultivariate{typeof(maniT)}(x, manifolds=maniT)
end

export ContinuousEuclid
struct ContinuousEuclid{N} <: InferenceVariable end

ContinuousEuclid(x::Int) = ContinuousEuclid{x}()

function Base.getproperty(::ContinuousEuclid{N}, f::Symbol) where N
  if f == :dims
    return N
  elseif f == :manifolds
    return ntuple(i -> :Euclid, N)
  end
end


"""
$(TYPEDEF)

Default prior on all dimensions of a variable node in the factor graph.  `Prior` is
not recommended when non-Euclidean dimensions are used in variables.
"""
struct Prior{T <: SamplableBelief} <: AbstractPrior 
  Z::T
end
getSample(s::Prior, N::Int=1) = (reshape(rand(s.Z,N),:,N), )


# TODO maybe replace X with a type.
function (s::Prior{<:ParametricTypes})(X1::AbstractVector{T};
                    userdata::Union{Nothing,FactorMetadata}=nothing) where T <: Real

  if isa(s.Z, Normal)
    meas = s.Z.μ
    σ = s.Z.σ
    #TODO confirm signs
    res = meas - X1[1]
    return (res./σ) .^ 2

  elseif isa(s.Z, MvNormal)
    meas = mean(s.Z)
    iΣ = invcov(s.Z)
    #TODO confirm math : Σ^(1/2)*X
    res = meas .- X1
    return res' * iΣ * res # + 2*log(1/(  sqrt(det(Σ)*(2pi)^k) )) ## cancel ×1/2 in calling function ## k = dim(μ)
  else
    #this should not happen
    @error("$s not suported, please use non-parametric")
  end
end

"""
$(TYPEDEF)

Message prior on all dimensions of a variable node in the factor graph.

Notes
- Only temporary existance during CSM operations.
"""
struct MsgPrior{T <: SamplableBelief} <: AbstractPrior
  Z::T
  inferdim::Float64
  MsgPrior{T}() where {T} = new{T}()
  MsgPrior{T}(z::T, infd::R) where {T <: SamplableBelief, R <: Real} = new{T}(z, infd)
end
function MsgPrior(z::T, infd::R) where {T <: SamplableBelief, R <: Real}
    MsgPrior{T}(z, infd)
end
getSample(s::MsgPrior, N::Int=1) = (reshape(rand(s.Z,N),:,N), )

function (s::MsgPrior{<:ParametricTypes})(X1::AbstractVector{T};
                       userdata::Union{Nothing,FactorMetadata}=nothing) where T<:Real

  if isa(s.Z, Normal)
    meas = s.Z.μ
    σ = s.Z.σ
    #TODO confirm signs
    res = meas - X1[1]
    return (res./σ) .^ 2

  elseif isa(s.Z, MvNormal)
    meas = mean(s.Z)
    iΣ = invcov(s.Z)
    #TODO confirm math : Σ^(1/2)*X
    res = meas .- X1
    return res' * iΣ * res

  else
    #this should not happen
    @error("$s not suported, please use non-parametric")
  end                    #
end

struct PackedMsgPrior <: PackedInferenceType
  Z::String
  inferdim::Float64
  PackedMsgPrior() = new()
  PackedMsgPrior(z::S, infd::R) where {S <: AbstractString, R <: Real} = new(string(z), infd)
end

function convert(::Type{PackedMsgPrior}, d::MsgPrior)
  PackedMsgPrior(string(d.Z), d.inferdim)
end
function convert(::Type{MsgPrior}, d::PackedMsgPrior)
  MsgPrior(extractdistribution(d.Z), d.inferdim)
end


"""
Converter: Prior -> PackedPrior::Dict{String, Any}
"""
function convert(::Type{Dict{String, Any}}, prior::IncrementalInference.Prior)
    z = convert(Type{Dict{String, Any}}, prior.Z)
    return Packed_Factor([z], "Prior")
end

"""
Converter: PackedPrior::Dict{String, Any} -> Prior
"""
function convert(::Type{Prior}, prior::Dict{String, Any})
    # Genericize to any packed type next.
    z = prior["measurement"][1]
    z = convert(_evalType(z["distType"]), z)
    return Prior(z)
end

"""
$(TYPEDEF)

Partial prior belief (absolute data) on any variable, given `<:SamplableBelief` and which dimensions of the intended variable.
"""
struct PartialPrior{T <: SamplableBelief,P <: Tuple} <: AbstractPrior
  Z::T
  partial::P
end
getSample(s::PartialPrior, N::Int=1) = (reshape(rand(s.Z,N),:,N), )



"""
$(TYPEDEF)

Default linear offset between two scalar variables.
"""
struct LinearConditional{N, T <: SamplableBelief} <: AbstractRelativeFactor
  Z::T
end
function LinearConditional{N}() where N
  newval = MvNormal(zeros(N), diagm(ones(N)))
  LinearConditional{N,typeof(newval)}(newval)
end
LinearConditional() = LinearConditional{1}()
LinearConditional(nm::Distributions.ContinuousUnivariateDistribution) = LinearConditional{1, typeof(nm)}(nm)
LinearConditional(nm::BallTreeDensity) = LinearConditional{Ndim(nm), typeof(nm)}(nm)


getSample(s::LinearConditional, N::Int=1) = (reshape(rand(s.Z,N),:,N), )
function (s::LinearConditional)(res::AbstractArray{<:Real},
                                userdata::FactorMetadata,
                                idx::Int,
                                meas::Tuple{<:AbstractArray{<:Real, 2}},
                                X1::AbstractArray{<:Real,2},
                                X2::AbstractArray{<:Real,2}  )
  #
  #   @show size(res), size(meas[1]), size(X1), size(X2)
  res[:] = meas[1][:,idx] - (X2[:,idx] - X1[:,idx])
  nothing
end

# parametric specific functor
function (s::LinearConditional{<:ParametricTypes})(
                                X1::AbstractVector{<:Real},
                                X2::AbstractVector{<:Real};
                                userdata::Union{Nothing,FactorMetadata}=nothing )
                                #can I change userdata to a keyword arg
  #
  # FIXME, replace if with dispatch
  if isa(s.Z, Normal)
    meas = mean(s.Z)
    σ = std(s.Z)
    # res = similar(X2)
    res = meas - (X2[1] - X1[1])
    return (res/σ) .^ 2

  elseif isa(s.Z, MvNormal)
    meas = mean(s.Z)
    iΣ = invcov(s.Z)
    #TODO confirm math : Σ^(1/2)*X
    res = meas .- (X2 .- X1)
    return res' * iΣ * res

  else
    #this should not happen
    @error("$s not suported, please use non-parametric")
  end
end


"""
$(TYPEDEF)

Define a categorical mixture of prior beliefs on a variable.
"""
mutable struct MixturePrior{T <: SamplableBelief} <: AbstractPrior
  Z::Vector{T}
  C::Distributions.Categorical
  #derived values
  labels::Vector{Int}
  dims::Ref{Int}
  smpls::Array{Float64,2}
  MixturePrior{T}() where T  = new{T}()
  MixturePrior{T}(z::Vector{T}, c::DiscreteNonParametric) where {T <: SamplableBelief} = new{T}(z, c, zeros(Int,0),0,zeros(1,0))
  MixturePrior{T}(z::Vector{T}, p::Vector{Float64}) where {T <: SamplableBelief} = new{T}(z, Distributions.Categorical(p), zeros(Int,0),0,zeros(1,0))
end

function MixturePrior(z::Vector{T}, c::Union{<:Distributions.DiscreteNonParametric, Vector{Float64}})::MixturePrior{T} where {T <: SamplableBelief}
  MixturePrior{T}(z, c)
end

function getSample(s::MixturePrior, N::Int=1)
  #out memory should be right size first
  (length(s.labels) != N) && resize!(s.labels, N)
  (s.dims[] != size(s.smpls,1)) && ( s.dims[] = size( rand(s.Z[1],1), 1) )
  (size(s.smpls,2) != N) && ( s.smpls = Array{Float64,2}(undef,s.dims[],N) )
  s.labels .= rand(s.C, N)
  for i in 1:N
    s.smpls[:,i] .= rand(s.Z[s.labels[i]],1)
  end
  (s.smpls, s.labels)
end


"""
$(TYPEDEF)

Define a categorical mixture of (relative) likelihood beliefs between any two variables.
"""
struct MixtureLinearConditional{T} <: AbstractRelativeFactor
  Z::Vector{T}
  C::Distributions.Categorical
  MixtureLinearConditional{T}() where T  = new{T}()
  MixtureLinearConditional{T}(z::Vector{T}, c::Distributions.Categorical) where {T <: SamplableBelief} = new{T}(z, c)
  MixtureLinearConditional{T}(z::Vector{T}, p::Vector{Float64}) where {T <: SamplableBelief} = MixtureLinearConditional{T}(z, Distributions.Categorical(p))
end
MixtureLinearConditional(z::Vector{T}, c::Union{Distributions.Categorical, Vector{Float64}}) where {T <: SamplableBelief} = MixtureLinearConditional{T}(z, c)

getSample(s::MixtureLinearConditional, N::Int=1) = (reshape.(rand.(s.Z, N),1,:)..., rand(s.C, N))
function (s::MixtureLinearConditional)(res::AbstractArray{<:Real},
                               userdata::FactorMetadata,
                               idx::Int,
                               meas::Tuple,
                               X1::AbstractArray{<:Real,2},
                               X2::AbstractArray{<:Real,2}  )
  #
  res[1] = meas[meas[end][idx]][idx] - (X2[1,idx] - X1[1,idx])
  nothing
end



## packed types are still developed by hand.  Future versions would likely use a @packable macro to write Protobuf safe versions of factors

"""
$(TYPEDEF)

Serialization type for Prior.
"""
mutable struct PackedPrior <: PackedInferenceType
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


"""
$(TYPEDEF)

Serialization type for `PartialPrior`.
"""
mutable struct PackedPartialPrior <: PackedInferenceType
  Z::String
  partials::Vector{Int}
  PackedPartialPrior() = new()
  PackedPartialPrior(z::AS) where {AS <: AbstractString} = new(z)
end
function convert(::Type{PackedPartialPrior}, d::PartialPrior)
  PackedPartialPrior(string(d.Z), [d.partial...;])
end
function convert(::Type{PartialPrior}, d::PackedPartialPrior)
  PartialPrior(extractdistribution(d.Z),(d.partials...,))
end


"""
$(TYPEDEF)

Serialization type for `LinearConditional` binary factor.
"""
mutable struct PackedLinearConditional <: PackedInferenceType
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


"""
$(TYPEDEF)

Serialization type for `MixtureLinearConditional`.
"""
mutable struct PackedMixtureLinearConditional <: PackedInferenceType
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



"""
$(TYPEDEF)

Serialization type for `MixedPrior`.
"""
mutable struct PackedMixturePrior <: PackedInferenceType
  strs::Vector{String}
  cat::String
  PackedMixturePrior() = new()
  PackedMixturePrior(z::Vector{<:AbstractString}, cstr::AS) where {AS <: AbstractString} = new(z, cstr)
end
function convert(::Type{PackedMixturePrior}, d::MixturePrior)
  PackedMixturePrior(string.(d.Z), string(d.C))
end
function convert(::Type{MixturePrior}, d::PackedMixturePrior)
  MixturePrior(extractdistribution.(d.strs), extractdistribution(d.cat))
end
