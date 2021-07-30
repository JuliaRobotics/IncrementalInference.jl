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
    # return log(M, q, q̂)
    return vee(M, q, log(M, q, q̂))
    # return distance(M, q, q̂)
end

#::MeasurementOnTangent
function distanceTangent2Point(M::AbstractGroupManifold, X, p, q)
    q̂ = compose(M, p, exp(M, identity(M, p), X)) #for groups
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
abstract type AbstractManifoldMinimize <: AbstractRelative end

export ManifoldFactor
# DEV NOTES
# For now, `Z` is on the tangent space in coordinates at the point used in the factor.
# For groups just the lie algebra
# As transition it will be easier this way, we can reevaluate
struct ManifoldFactor{M <: AbstractManifold, T <: SamplableBelief} <: AbstractManifoldMinimize#AbstractFactor
    M::M
    Z::T
end

DFG.getManifold(f::ManifoldFactor) = f.M

# function getSample(cf::ManifoldFactor, N::Int=1)
function getSample(cf::CalcFactor{<:ManifoldFactor}, N::Int=1)
    #TODO @assert dim == cf.factor.Z's dimension
    #TODO investigate use of SVector if small dims
    ret = [rand(cf.factor.Z) for _ in 1:N]

    #TODO tangent or not?
    # tangent for now to fit with rest
    (ret, )
end

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
struct ManifoldPrior{M <: AbstractManifold, T <: SamplableBelief, P} <: AbstractPrior
    M::M 
    p::P #NOTE This is a fixed point from where the measurement `Z` is made in coordinates on tangent TpM
    Z::T
end

DFG.getManifold(f::ManifoldPrior) = f.M

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
    
    # X = get_vector.(Ref(M), Ref(p), Xc, Ref(DefaultOrthogonalBasis()))
    X = hat.(Ref(M), Ref(p), Xc)
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

if false
using IncrementalInference
using Manifolds
using LinearAlgebra
using StaticArrays

f = ManifoldFactor(SpecialOrthogonal(3), MvNormal([0.1, 0.02, 0.01]))
s = getSample(f,10)[1]
s[1]

f = ManifoldFactor(SpecialEuclidean(2), MvNormal([0.1, 0.2, 0.01]))
s = getSample(f,10)[1]
s[1]


f = ManifoldPrior(SpecialOrthogonal(2), SA[1.0 0; 0 1], MvNormal([0.1]))
meas = getSample(f,10)[1]
meas[1]
f.(meas, Ref(SA[1.0 0; 0 1]))

f = ManifoldPrior(SpecialOrthogonal(3), SA[1.0 0 0; 0 1 0; 0 0 1], MvNormal([0.1, 0.02, 0.01]))
s = getSample(f,10)[1]
s[1]

f = ManifoldPrior(SpecialEuclidean(2), ProductRepr(SA[0,0], SA[1.0 0; 0 1]), MvNormal([0.1, 0.2, 0.01]))
s = getSample(f,10)[1]
s[1]



end
