

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






##

"""
    $SIGNATURES

First hacky version to return which factor type to use between two variables of types T1 and T2.
"""
selectFactorType(T1::Type{ContinuousScalar}, T2::Type{ContinuousScalar}) = LinearRelative
selectFactorType(T1::Type{<:ContinuousEuclid{N}}, T2::Type{<:ContinuousEuclid{N}}) where N = LinearRelative{N}
selectFactorType(T1::InferenceVariable, T2::InferenceVariable) = selectFactorType(typeof(T1), typeof(T2))
selectFactorType(dfg::AbstractDFG, s1::Symbol, s2::Symbol) = selectFactorType( getVariableType(dfg, s1), getVariableType(dfg, s2) )

