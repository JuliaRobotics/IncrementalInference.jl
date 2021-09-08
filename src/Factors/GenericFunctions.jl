# TODO under development - experimenting with type to work with manifolds
## ======================================================================================
## Generic manifold cost functions
## ======================================================================================
"""
    $SIGNATURES
Generic function that can be used in binary factors to calculate distance between points on Lie Groups with measurements.
"""
function distancePoint2Point(M::AbstractGroupManifold, m, p, q)
    q̂ = Manifolds.compose(M, p, m)
    # return log(M, q, q̂)
    return vee(M, q, log(M, q, q̂))
    # return distance(M, q, q̂)
end

#::MeasurementOnTangent
function distanceTangent2Point(M::AbstractGroupManifold, X, p, q)
    q̂ = Manifolds.compose(M, p, exp(M, identity_element(M, p), X)) #for groups
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
abstract type AbstractManifoldMinimize <: AbstractRelative end

export ManifoldFactor
# DEV NOTES
# For now, `Z` is on the tangent space in coordinates at the point used in the factor.
# For groups just the lie algebra
# As transition it will be easier this way, we can reevaluate
struct ManifoldFactor{M <: AbstractManifold, T <: SamplableBelief} <: AbstractManifoldMinimize #AbstractFactor
    M::M
    Z::T
end

DFG.getManifold(f::ManifoldFactor) = f.M

function getSample(cf::CalcFactor{<:ManifoldFactor{M,Z}}) where {M,Z}
    #TODO @assert dim == cf.factor.Z's dimension
    #TODO investigate use of SVector if small dims
    if M isa ManifoldKernelDensity
        ret = sample(cf.factor.Z.belief)[1]
    else
        ret = rand(cf.factor.Z)
    end
    #return coordinates as we do not know the point here #TODO separate Lie group
    return ret
end

# function (cf::CalcFactor{<:ManifoldFactor{<:AbstractGroupManifold}})(Xc, p, q)
function (cf::CalcFactor{<:ManifoldFactor})(Xc, p, q)
# function (cf::ManifoldFactor)(X, p, q)
    M = cf.factor.M
    # M = cf.M
    X = hat(M, p, Xc)
    return distanceTangent2Point(M, X, p, q)
end

## ======================================================================================
## ManifoldPrior
## ======================================================================================
export ManifoldPrior
# `p` is a point on manifold `M`
# `Z` is a measurement at the tangent space of `p` on manifold `M` 
struct ManifoldPrior{M <: AbstractManifold, T <: SamplableBelief, P, B <: AbstractBasis} <: AbstractPrior
    M::M 
    p::P #NOTE This is a fixed point from where the measurement `Z` is made in coordinates on tangent TpM
    Z::T
    basis::B
    retract_method::AbstractRetractionMethod
end

ManifoldPrior(M::AbstractGroupManifold, p, Z) = ManifoldPrior(M, p, Z, ManifoldsBase.VeeOrthogonalBasis(), ExponentialRetraction())

DFG.getManifold(f::ManifoldPrior) = f.M

#TODO
# function ManifoldPrior(M::AbstractGroupManifold, Z::SamplableBelief)
#     # p = identity_element(M, #TOOD)
#     # similar to getPointIdentity(M)
#     return ManifoldPrior(M, Z, p)
# end

# ManifoldPrior{M}(Z::SamplableBelief, p) where M = ManifoldPrior{M, typeof(Z), typeof(p)}(Z, p)

function getSample(cf::CalcFactor{<:ManifoldPrior})
    Z = cf.factor.Z
    p = cf.factor.p
    M = cf.factor.M
    basis = cf.factor.basis
    retract_method = cf.factor.retract_method

    point = samplePoint(M, Z, p, basis, retract_method)
    
    return point
end

#TODO investigate SVector if small dims, this is slower
# dim = manifold_dimension(M)
# Xc = [SVector{dim}(rand(Z)) for _ in 1:N]

function (cf::CalcFactor{<:ManifoldPrior})(m, p)
    M = cf.factor.M
    # return log(M, p, m)
    return vee(M,p,log(M, p, m))
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
