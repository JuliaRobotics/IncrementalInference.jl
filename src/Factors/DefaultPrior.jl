
export ManifoldPrior



"""
$(TYPEDEF)

Default prior on all dimensions of a variable node in the factor graph.  `Prior` is
not recommended when non-Euclidean dimensions are used in variables.
"""
struct Prior{T <: SamplableBelief} <: AbstractPrior 
  Z::T
end
Prior(::UniformScaling) = Prior(Normal())

function getSample(cf::CalcFactor{<:Prior}, N::Int=1)
  meas = Vector{Vector{Float64}}(undef, N)
  for i in 1:N
    meas[i] = rand(cf.factor.Z,1)[:]
  end
  (meas,)
end

# basic default
(s::CalcFactor{<:Prior})(z, x1) = z .- x1

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
  PackedPrior(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{Prior}, d::PackedPrior)
  Prior(convert(SamplableBelief, d.Z))
end



"""
Converter: Prior -> PackedPrior::Dict{String, Any}

FIXME see DFG #590 for consolidation with Serialization and Marshaling
"""
function convert(::Type{Dict{String, Any}}, prior::IncrementalInference.Prior)
    z = convert(Type{Dict{String, Any}}, prior.Z)
    return Packed_Factor([z], "Prior")
end

"""
Converter: PackedPrior::Dict{String, Any} -> Prior

FIXME see DFG #590 for consolidation on Serialization and Marshaling
"""
function convert(::Type{<:Prior}, prior::Dict{String, Any})
    # Genericize to any packed type next.
    z = prior["measurement"][1]
    z = convert(_evalType(z["distType"]), z)
    return Prior(z)
end



## ======================================================================================
## ManifoldPrior
## ======================================================================================
# `p` is a point on manifold `M`
# `Z` is a measurement at the tangent space of `p` on manifold `M` 
# PLEASE NOTE, ManifoldPrior CANNOT BE SERIALIZED, DO NOT USE!
struct ManifoldPrior{M <: AbstractManifold, T <: SamplableBelief, P} <: AbstractPrior
  M::M 
  Z::T
  p::P
end

#TODO
# function ManifoldPrior(M::AbstractGroupManifold, Z::SamplableBelief)
#     # p = identity(M, #TOOD)
#     # similar to getPointIdentity(M)
#     return ManifoldPrior(M, Z, p)
# end

# ManifoldPrior{M}(Z::SamplableBelief, p) where M = ManifoldPrior{M, typeof(Z), typeof(p)}(Z, p)

# function getSample(cf::ManifoldPrior, N::Int=1)
function getSample(cf::CalcFactor{<:ManifoldPrior}, N::Int=1)
  Z = cf.factor.Z
  p = cf.factor.p
  M = cf.factor.M
  # Z = cf.Z
  # p = cf.p
  # M = cf.M
  
  Xc = [rand(Z) for _ in 1:N]
  
  X = get_vector.(Ref(M), Ref(p), Xc, Ref(DefaultOrthogonalBasis()))
  points = exp.(Ref(M), Ref(p), X)

  return (points, )
end

#TODO investigate SVector if small dims, this is slower
# dim = manifold_dimension(M)
# Xc = [SVector{dim}(rand(Z)) for _ in 1:N]

# function (cf::ManifoldPrior)(m, p)
function (cf::CalcFactor{<:ManifoldPrior})(m, p)
  M = cf.factor.M
  # M = cf.M
  return distancePrior(M, m, p)
end
