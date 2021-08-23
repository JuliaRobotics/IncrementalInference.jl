
export ContinuousScalar
export ContinuousEuclid
export Circular, Circle


Base.convert(::Type{<:Tuple}, ::InstanceType{Manifolds.Euclidean{Tuple{N}, ℝ}} ) where N = tuple([:Euclid for i in 1:N]...)
Base.convert(::Type{<:Tuple}, ::InstanceType{Manifolds.Circle{ℝ}})  = (:Circular,)




## Euclid 1

"""
$(TYPEDEF)

Most basic continuous scalar variable in a `::DFG.AbstractDFG` object.

DevNotes
- TODO Consolidate with ContinuousEuclid{1}
"""
@defVariable ContinuousScalar TranslationGroup(1) [0.0;]



"""
    ContinuousEuclid{N}
Continuous Euclidean variable of dimension `N`.
"""
struct ContinuousEuclid{N} <: InferenceVariable end

ContinuousEuclid(x::Int) = ContinuousEuclid{x}()

# not sure if these overloads are necessary since DFG 775?
DFG.getManifold(::InstanceType{ContinuousEuclid{N}}) where N = TranslationGroup(N)                             
DFG.getDimension(val::InstanceType{ContinuousEuclid{N}}) where N = manifold_dimension(getManifold(val))


DFG.getPointType(::Type{ContinuousEuclid{N}}) where N = Vector{Float64}
DFG.getPointIdentity(M_::Type{ContinuousEuclid{N}}) where N = zeros(N) # identity_element(getManifold(M_), zeros(N)) 


Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::InstanceType{ContinuousEuclid{N}}) where N = TranslationGroup(N)


## Circular


"""
$(TYPEDEF)

Circular is a `Manifolds.Circle{ℝ}` mechanization of one rotation, with `theta in [-pi,pi)`.
"""
@defVariable Circular Circle() [0.0;]




#