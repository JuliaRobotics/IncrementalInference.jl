
export ContinuousScalar
export ContinuousEuclid
export Circular, Circle


Base.convert(::Type{<:Tuple}, ::InstanceType{Manifolds.Euclidean{Tuple{N}, ℝ}} ) where N = tuple([:Euclid for i in 1:N]...)
Base.convert(::Type{<:Tuple}, ::InstanceType{Manifolds.Circle{ℝ}})  = (:Circular,)

# Base.convert(::Type{<:Tuple}, mani::ManifoldsBase.Manifold) = getManifolds(mani)


## Euclid 1

"""
$(TYPEDEF)

Most basic continuous scalar variable in a `::DFG.AbstractDFG` object.

DevNotes
- TODO Consolidate with ContinuousEuclid{1}
"""
@defVariable ContinuousScalar Euclidean(1) #1 (:Euclid,)


## Euclid N


"""
    ContinuousEuclid{N}
Continuous Euclidean variable of dimension `N`.
"""
struct ContinuousEuclid{N} <: InferenceVariable end

ContinuousEuclid(x::Int) = ContinuousEuclid{x}()

getManifold(::Type{<:ContinuousEuclid{N}}) where N = Euclidean(N)                               # ntuple(i -> :Euclid, N)
getDimension(val::Type{<:ContinuousEuclid{N}}) where N = manifold_dimension(getManifold(val))   # N::Int
getManifolds(val::Type{<:ContinuousEuclid{N}}) where N = convert(Tuple, getManifold(val))       # ntuple(i -> :Euclid, N)

getManifold(::ContinuousEuclid{N}) where N = Euclidean(N)                               # ntuple(i -> :Euclid, N)
getDimension(val::ContinuousEuclid{N}) where N = manifold_dimension(getManifold(val))   # N::Int
getManifolds(val::ContinuousEuclid{N}) where N = convert(Tuple, getManifold(val))       # ntuple(i -> :Euclid, N)

Base.convert(::Type{<:ManifoldsBase.Manifold}, ::InstanceType{ContinuousEuclid{N}}) where N = Manifolds.Euclidean(N)
Base.convert(::Type{<:ManifoldsBase.Manifold}, ::InstanceType{ContinuousScalar})    = Manifolds.Euclidean(1)


## Circular


"""
$(TYPEDEF)

Circular is a `Manifolds.Circle{ℝ}` mechanization of one rotation, with `theta in [-pi,pi)`.
"""
@defVariable Circular Circle()


Base.convert(::Type{<:ManifoldsBase.Manifold}, ::InstanceType{Circular}) = Manifolds.Circle()



#