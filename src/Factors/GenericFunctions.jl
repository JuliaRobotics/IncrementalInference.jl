# TODO under development - experimenting with type to work with manifolds
## ======================================================================================
## Generic manifold cost functions
## ======================================================================================
"""
    $SIGNATURES
Generic function that can be used in binary factors to calculate distance between points on Lie Groups with measurements.
"""
function distancePoint2Point(M::AbstractGroupManifold, m, p, q)
    q̂ = compose(M, p, m)
    return distance(M, q, q̂)
end

#::MeasurementOnTangent
function distanceTangent2Point(M::AbstractGroupManifold, X, p, q)
    q̂ = compose(M, p, exp(M, identity(M, p), X)) #for groups
    return distance(M, q, q̂)
end

# ::MeasurementOnTangent
function distanceTangent2Point(M::AbstractManifold, X, p, q)
    q̂ = exp(M, p, X) 
    return distance(M, q, q̂)
end

"""
    $SIGNATURES
Generic function that can be used in prior factors to calculate distance on Lie Groups. 
"""
function distancePrior(M::AbstractManifold, meas, p)	
    return distance(M, meas, p)
end

## ======================================================================================
## Default Identities #TODO only development, replace with better idea
## ======================================================================================

default_identity(M) = error("No default identity element defined for $(typeof(M))")
default_identity(M::GroupManifold{ℝ, <:ProductManifold}) = error("No default identity element defined for $(typeof(M))")
function default_identity(::SpecialEuclidean{N}) where N
    T = Float64
    t = zeros(SVector{N, T})
    R = SMatrix{N,N,T}(one(T)I)
    return ProductRepr(t, R)
end
function default_identity(M::GroupManifold{ℝ, <:AbstractManifold})
    T = Float64
    s = representation_size(M)
    return identity(M, zeros(SArray{Tuple{s...},T}))
end

# function default_identity(::SpecialOrthogonal{N}) where N
#     T = Float64
#     return SMatrix{N,N,T}(one(T)I)
# end

## ======================================================================================
## ManifoldFactor
## ======================================================================================


export ManifoldFactor
# DEV NOTES
# For now, `Z` is on the tangent space in coordinates at the point used in the factor.
# For groups just the lie algebra
# As transition it will be easier this way, we can reevaluate
struct ManifoldFactor{M <: AbstractManifold, T <: SamplableBelief} <: AbstractManifoldMinimize#AbstractFactor
    M::M
    Z::T
end

# function getSample(cf::ManifoldFactor, N::Int=1)
function getSample(cf::CalcFactor{<:ManifoldFactor}, N::Int=1)
    #TODO @assert dim == cf.factor.Z's dimension
    #TODO investigate use of SVector if small dims
    ret = [rand(cf.factor.Z) for _ in 1:N]

    #TODO tangent or not?
    # tangent for now to fit with rest
    (ret, )
end

function (cf::CalcFactor{<:ManifoldFactor})(X, p, q)
# function (cf::ManifoldFactor)(X, p, q)
    M = cf.factor.M
    # M = cf.M
    return distanceTangent2Point(M, X, p, q)
end


