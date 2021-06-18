
export ContinuousScalar
export ContinuousEuclid
export Circular, Circle


Base.convert(::Type{<:Tuple}, ::InstanceType{Manifolds.Euclidean{Tuple{N}, ℝ}} ) where N = tuple([:Euclid for i in 1:N]...)
Base.convert(::Type{<:Tuple}, ::InstanceType{Manifolds.Circle{ℝ}})  = (:Circular,)

# Base.convert(::Type{<:Tuple}, mani::MB.AbstractManifold) = getManifolds(mani)


## Euclid 1

"""
$(TYPEDEF)

Most basic continuous scalar variable in a `::DFG.AbstractDFG` object.

DevNotes
- TODO Consolidate with ContinuousEuclid{1}
"""
@defVariable ContinuousScalar Euclidean(1) SVector{1}

# Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::InstanceType{ContinuousScalar})    = Manifolds.Euclidean(1)


## Euclid N


"""
    ContinuousEuclid{N}
Continuous Euclidean variable of dimension `N`.
"""
struct ContinuousEuclid{N} <: InferenceVariable end

ContinuousEuclid(x::Int) = ContinuousEuclid{x}()

getManifold(::Type{<:ContinuousEuclid{N}}) where N = Euclidean(N)
getDimension(val::Type{<:ContinuousEuclid{N}}) where N = manifold_dimension(getManifold(val))
# getManifolds(val::Type{<:ContinuousEuclid{N}}) where N = convert(Tuple, getManifold(val))

getManifold(::ContinuousEuclid{N}) where N = Euclidean(N)                               
getDimension(val::ContinuousEuclid{N}) where N = manifold_dimension(getManifold(val))
# getManifolds(val::ContinuousEuclid{N}) where N = convert(Tuple, getManifold(val))

Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::InstanceType{ContinuousEuclid{N}}) where N = Manifolds.Euclidean(N)


## Circular


"""
$(TYPEDEF)

Circular is a `Manifolds.Circle{ℝ}` mechanization of one rotation, with `theta in [-pi,pi)`.
"""
@defVariable Circular Circle() SVector{1}


# Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::InstanceType{Circular}) = Manifolds.Circle()



#