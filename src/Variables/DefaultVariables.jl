
export ContinuousEuclid


Base.convert(::Type{<:Tuple}, mani::ManifoldsBase.Manifold) = getManifolds(mani)


"""
$(TYPEDEF)

Most basic continuous scalar variable in a `::DFG.AbstractDFG` object.

DevNotes
- TODO Consolidate with ContinuousEuclid{1}
"""
@defVariable ContinuousScalar 1 (:Euclid,)

"""
$(TYPEDEF)

Continuous variable of dimension `.dims` on manifold `.manifolds`.
"""
struct ContinuousMultivariate{T1 <: Tuple} <: InferenceVariable
  dims::Int
  manifolds::T1
end

function ContinuousMultivariate(x::Int;
                                manifolds::T1=convert(Tuple, Euclidean(1))  )  where {T1 <: Tuple}
  #
  maniT = length(manifolds) < x ? ([manifolds[1] for i in 1:x]...,) : manifolds
  ContinuousMultivariate{typeof(maniT)}(x, maniT)
end


"""
    ContinuousEuclid{N}
Continuous Euclidean variable of dimension `N`.
"""
struct ContinuousEuclid{N} <: InferenceVariable end

ContinuousEuclid(x::Int) = ContinuousEuclid{x}()

getDimension(::ContinuousEuclid{N}) where N = N::Int
getManifolds(::ContinuousEuclid{N}) where N = ntuple(i -> :Euclid, N)


Base.convert(::Type{<:ManifoldsBase.Manifold}, ::InstanceType{ContinuousEuclid{N}}) where N = Manifolds.Euclidean(N)
Base.convert(::Type{<:ManifoldsBase.Manifold}, ::InstanceType{ContinuousScalar})    = Manifolds.Euclidean(1)

