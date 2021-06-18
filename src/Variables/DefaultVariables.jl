
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
@defVariable ContinuousScalar Euclidean(1) Vector{Float64}

# Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::InstanceType{ContinuousScalar})    = Manifolds.Euclidean(1)


## Euclid N


"""
    ContinuousEuclid{N}
Continuous Euclidean variable of dimension `N`.
"""
struct ContinuousEuclid{N} <: InferenceVariable end

ContinuousEuclid(x::Int) = ContinuousEuclid{x}()
DFG.getPointType(::Type{ContinuousEuclid{N}}) where N = Vector{Float64}

# not sure if these overloads are necessary since DFG 775?
DFG.getManifold(::Type{<:ContinuousEuclid{N}}) where N = Euclidean(N)
DFG.getDimension(val::Type{<:ContinuousEuclid{N}}) where N = manifold_dimension(getManifold(val))
DFG.getManifold(::ContinuousEuclid{N}) where N = Euclidean(N)                               
DFG.getDimension(val::ContinuousEuclid{N}) where N = manifold_dimension(getManifold(val))


Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::InstanceType{ContinuousEuclid{N}}) where N = Manifolds.Euclidean(N)


## Circular


"""
$(TYPEDEF)

Circular is a `Manifolds.Circle{ℝ}` mechanization of one rotation, with `theta in [-pi,pi)`.
"""
@defVariable Circular Circle() Vector{Float64}


# Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::InstanceType{Circular}) = Manifolds.Circle()



#