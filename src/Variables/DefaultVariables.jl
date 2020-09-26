
export ContinuousEuclid


"""
$(TYPEDEF)

Most basic continuous scalar variable in a `::DFG.AbstractDFG` object.
"""
struct ContinuousScalar <: InferenceVariable
  function ContinuousScalar(;manifolds=nothing)
    manifolds !== nothing &&
    Base.depwarn("ContinuousScalar keyword argument manifolds is deprecated.", :ContinuousScalar)
    return new()
  end
end
getDimension(::ContinuousScalar) = 1
getManifolds(::ContinuousScalar) = (:Euclid,)

"""
$(TYPEDEF)

Continuous variable of dimension `.dims` on manifold `.manifolds`.
"""
struct ContinuousMultivariate{T1 <: Tuple} <: InferenceVariable
  dims::Int
  manifolds::T1
  ContinuousMultivariate{T}() where {T} = new()
  ContinuousMultivariate{T}(x::Int; manifolds::T=(:Euclid,)) where {T <: Tuple} = new(x, manifolds)
end

function ContinuousMultivariate(x::Int;
                                manifolds::T1=(:Euclid,)  )  where {T1 <: Tuple}
  #
  maniT = length(manifolds) < x ? ([manifolds[1] for i in 1:x]...,) : manifolds
  ContinuousMultivariate{typeof(maniT)}(x, manifolds=maniT)
end


"""
    ContinuousEuclid{N}
Continuous Euclidean variable of dimension `N`.
"""
struct ContinuousEuclid{N} <: InferenceVariable end

ContinuousEuclid(x::Int) = ContinuousEuclid{x}()

getDimension(::ContinuousEuclid{N}) where N = N::Int
getManifolds(::ContinuousEuclid{N}) where N = ntuple(i -> :Euclid, N)
