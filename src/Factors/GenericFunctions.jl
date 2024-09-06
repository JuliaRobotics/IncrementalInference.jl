# TODO under development - experimenting with type to work with manifolds

## ======================================================================================
## Possible type piracy, but also standardizing to common API across repos
## ======================================================================================

DFG.getDimension(::Distributions.Uniform) = 1
DFG.getDimension(::Normal) = 1
DFG.getDimension(Z::MvNormal) = Z |> cov |> diag |> length

function DFG.getDimension(Z::FluxModelsDistribution)
  return if length(Z.outputDim) == 1
    Z.outputDim[1]
  else
    error(
    "can only do single index tensor at this time, please open an issue with Caesar.jl",
  )
  end
end
DFG.getDimension(Z::ManifoldKernelDensity) = getManifold(Z) |> getDimension
# TODO deprecate
DFG.getDimension(Z::BallTreeDensity) = Ndim(Z)

## ======================================================================================
## Generic manifold cost functions
## ======================================================================================
"""
    $SIGNATURES
Generic function that can be used in binary factors to calculate distance between points on Lie Groups with measurements.
"""
function distancePoint2Point(M::SemidirectProductGroup, m, p, q)
  q̂ = Manifolds.compose(M, p, m)
  # return log(M, q, q̂)
  return vee(M, q, log(M, q, q̂))
  # return distance(M, q, q̂)
end

#::MeasurementOnTangent
function distanceTangent2Point(M::SemidirectProductGroup, X, p, q)
  q̂ = Manifolds.compose(M, p, exp(M, getPointIdentity(M), X)) #for groups
  # return log(M, q, q̂)
  return vee(M, q, log(M, q, q̂))
  # return distance(M, q, q̂)
end

# ::MeasurementOnTangent
function distanceTangent2Point(M::AbstractManifold, X, p, q)
  q̂ = exp(M, p, X)
  # return log(M, q, q̂)
  return vee(M, q, log(M, q, q̂))
  # return distance(M, q, q̂)
end

"""
    $SIGNATURES
Generic function that can be used in prior factors to calculate distance on Lie Groups. 
"""
function distancePrior(M::AbstractManifold, meas, p)
  return log(M, p, meas)
  # return distance(M, p, meas)
end

## ======================================================================================
## ManifoldFactor
## ======================================================================================

export ManifoldFactor
# DEV NOTES
# For now, `Z` is on the tangent space in coordinates at the point used in the factor.
# For groups just the lie algebra
# As transition it will be easier this way, we can reevaluate
struct ManifoldFactor{
  M <: AbstractManifold, 
  T <: SamplableBelief
} <: AbstractManifoldMinimize
  M::M
  Z::T
end

DFG.getManifold(f::ManifoldFactor) = f.M

function getSample(cf::CalcFactor{<:ManifoldFactor{M, Z}}) where {M, Z}
  #TODO @assert dim == cf.factor.Z's dimension
  #TODO investigate use of SVector if small dims
  # if M isa ManifoldKernelDensity
  #   ret = sample(cf.factor.Z.belief)[1]
  # else
  #   ret = rand(cf.factor.Z)
  # end

  # ASSUME this function is only used for RelativeFactors which must use measurements as tangents
  ret = sampleTangent(cf.factor.M, cf.factor.Z)
  #return coordinates as we do not know the point here #TODO separate Lie group
  return ret
end

# function (cf::CalcFactor{<:ManifoldFactor{<:AbstractDecoratorManifold}})(Xc, p, q)
function (cf::CalcFactor{<:ManifoldFactor})(X, p, q)
  return distanceTangent2Point(cf.factor.M, X, p, q)
end


## ======================================================================================
## adjoint factor - adjoint action applied to the measurement
## ======================================================================================


# Adjoints defined in ApproxManifoldProducts
struct AdFactor{F <: AbstractManifoldMinimize} <: AbstractManifoldMinimize
  factor::F
end

function (cf::CalcFactor{<:AdFactor})(Xϵ, p, q)
  # M = getManifold(cf.factor)
  # p,q ∈ M
  # Xϵ ∈ TϵM
  # ϵ = identity_element(M)
  # transform measurement from TϵM to TpM (global to local coordinates)
  # Adₚ⁻¹ = AdjointMatrix(M, p)⁻¹ = AdjointMatrix(M, p⁻¹)
  # Xp = Adₚ⁻¹ * Xϵᵛ
  # ad = Ad(M, inv(M, p))
  # Xp = Ad(M, inv(M, p), Xϵ)
  # Xp = adjoint_action(M, inv(M, p), Xϵ)
  #TODO is vector transport supposed to be the same?
  # Xp = vector_transport_to(M, ϵ, Xϵ, p)

  # Transform measurement covariance
  # ᵉΣₚ = Adₚ ᵖΣₚ Adₚᵀ
  #TODO test if transforming sqrt_iΣ is the same as Σ
  # Σ = ad * inv(cf.sqrt_iΣ^2) * ad'
  # sqrt_iΣ = convert(typeof(cf.sqrt_iΣ), sqrt(inv(Σ)))
  # sqrt_iΣ = convert(typeof(cf.sqrt_iΣ), ad * cf.sqrt_iΣ * ad')
  Xp = Xϵ

  child_cf = CalcFactorResidual(
    cf.faclbl,
    cf.factor.factor,
    cf.varOrder,
    cf.varOrderIdxs,
    cf.meas,
    cf.sqrt_iΣ,
    cf.cache,
  )
  return child_cf(Xp, p, q)
end

getMeasurementParametric(f::AdFactor) = getMeasurementParametric(f.factor)

getManifold(f::AdFactor) = getManifold(f.factor)
function getSample(cf::CalcFactor{<:AdFactor}) 
  M = getManifold(cf)
  return sampleTangent(M, cf.factor.factor.Z)
end

## ======================================================================================
## ManifoldPrior
## ======================================================================================
export ManifoldPrior, PackedManifoldPrior

# `p` is a point on manifold `M`
# `Z` is a measurement at the tangent space of `p` on manifold `M` 
struct ManifoldPrior{M <: AbstractManifold, T <: SamplableBelief, P, B <: AbstractBasis} <:
       AbstractPrior
  M::M
  p::P #NOTE This is a fixed point from where the measurement `Z` is made in coordinates on tangent TpM
  Z::T
  basis::B
  retract_method::AbstractRetractionMethod
end

function ManifoldPrior(M::AbstractDecoratorManifold, p, Z)
  return ManifoldPrior(M, p, Z, ManifoldsBase.VeeOrthogonalBasis(), ExponentialRetraction())
end

DFG.getManifold(f::ManifoldPrior) = f.M

#TODO
# function ManifoldPrior(M::AbstractDecoratorManifold, Z::SamplableBelief)
#     # p = identity_element(M, #TOOD)
#     # similar to getPointIdentity(M)
#     return ManifoldPrior(M, Z, p)
# end

# ManifoldPrior{M}(Z::SamplableBelief, p) where M = ManifoldPrior{M, typeof(Z), typeof(p)}(Z, p)

function getSample(cf::CalcFactor{<:ManifoldPrior})
  Z = cf.factor.Z
  p = cf.factor.p
  M = cf.manifold # .factor.M
  basis = cf.factor.basis
  retract_method = cf.factor.retract_method
  point = samplePoint(M, Z, p, basis, retract_method)

  return point
end

function getFactorMeasurementParametric(fac::ManifoldPrior)
  M = getManifold(fac)
  dims = manifold_dimension(M)
  meas = fac.p
  iΣ = convert(SMatrix{dims, dims}, invcov(fac.Z))
  meas, iΣ
end

#TODO investigate SVector if small dims, this is slower
# dim = manifold_dimension(M)
# Xc = [SVector{dim}(rand(Z)) for _ in 1:N]

function (cf::CalcFactor{<:ManifoldPrior})(m, p)
  M = cf.factor.M
  # return log(M, p, m)
  return vee(M, p, log(M, p, m))
  # return distancePrior(M, m, p)
end

# dist²_Σ = ⟨X, Σ⁻¹*X'⟩
function mahalanobus_distance2(M, p, q, inv_Σ)
  Xc = log(M, p, q)
  return mahalanobus_distance2(Xc, inv_Σ)
end

function mahalanobus_distance2(M, X, inv_Σ)
  #TODO look to replace with inner(MM, p, X, inv_Σ*X)
  # Xc = get_coordinates(M, p, X, DefaultOrthogonalBasis())
  Xc = vee(M, p, X)
  return Xc' * inv_Σ * Xc
end

Base.@kwdef mutable struct PackedManifoldPrior <: AbstractPackedFactor
  varType::String
  p::Vector{Float64}  #NOTE This is a fixed point from where the measurement `Z` likely stored as a coordinate
  Z::PackedSamplableBelief
end

function convert(
  ::Union{Type{<:AbstractPackedFactor}, Type{<:PackedManifoldPrior}},
  obj::ManifoldPrior,
)
  #

  varT = DFG.typeModuleName(getVariableType(obj.M))

  c = AMP.makeCoordsFromPoint(obj.M, obj.p)

  # TODO convert all distributions to JSON
  Zst = convert(PackedSamplableBelief, obj.Z) # String

  return PackedManifoldPrior(varT, c, Zst)
end

function convert(
  ::Union{Type{<:AbstractFactor}, Type{<:ManifoldPrior}},
  obj::PackedManifoldPrior,
)
  #

  # piggy back on serialization of InferenceVariable rather than try serialize anything Manifolds.jl
  M = DFG.getTypeFromSerializationModule(obj.varType) |> getManifold

  # TODO this is too excessive
  e0 = getPointIdentity(M)
  # u0 = getPointIdentity(obj.varType)
  p = AMP.makePointFromCoords(M, obj.p, e0) #, u0)

  Z = convert(SamplableBelief, obj.Z)

  return ManifoldPrior(M, p, Z)
end

## ======================================================================================
## Generic Manifold Partial Prior
## ======================================================================================

function samplePointPartial(
  M::AbstractDecoratorManifold,
  z::Distribution,
  partial::Vector{Int},
  p = getPointIdentity(M),
  retraction_method::AbstractRetractionMethod = ExponentialRetraction(),
)
  dim = manifold_dimension(M)
  Xc = zeros(dim)
  Xc[partial] .= rand(z)
  X = hat(M, p, Xc)
  return retract(M, p, X, retraction_method)
end

struct ManifoldPriorPartial{M <: AbstractManifold, T <: SamplableBelief, P <: Tuple} <:
       AbstractPrior
  M::M
  Z::T
  partial::P
end

DFG.getManifold(f::ManifoldPriorPartial) = f.M

function getSample(cf::CalcFactor{<:ManifoldPriorPartial})
  Z = cf.factor.Z
  M = getManifold(cf)
  partial = collect(cf.factor.partial)

  return (samplePointPartial(M, Z, partial),)
end
