
##==============================================================================
## LEGACY SUPPORT FOR ZMQ IN CAESAR
##==============================================================================

export listSolvekeys

export _evalType

# not sure if and where this is still being used
function _evalType(pt::String)::Type
  @error "_evalType has been deprecated, use DFG serialization methods instead."
  try
    getfield(Main, Symbol(pt))
  catch ex
    io = IOBuffer()
    showerror(io, ex, catch_backtrace())
    err = String(take!(io))
    error("_evalType: Unable to locate factor/distribution type '$pt' in main context (e.g. do a using on all your factor libraries). Please check that this factor type is loaded into main. Stack trace = $err")
  end
end


##==============================================================================
## Deprecate code below before v0.25
##==============================================================================

# """
#     $(SIGNATURES)

# Multiply various full and partial dimension proposal densities.

# DevNotes
# - FIXME consolidate partial and full product AMP API, relates to #1010
# - TODO better consolidate with full dimension product
# - TODO -- reuse memory rather than rand here
# """
# function prodmultiplefullpartials(dens::Vector{BallTreeDensity},
#                                   partials::Dict{Any, Vector{BallTreeDensity}},
#                                   Ndims::Int,
#                                   N::Int,
#                                   manis::Tuple;
#                                   useExisting::Bool=false )
#   #
#   # calculate products over all dimensions, legacy proposals held in `dens` vector
#   pGM = AMP.manifoldProduct(dens, manis, Niter=1) |> getPoints

#   _partialProducts!(pGM, partials, manis, useExisting=useExisting)

#   return pGM
# end

# function _setCCWDecisionDimsConv!(ccwl::Union{CommonConvWrapper{F},
#                                               CommonConvWrapper{Mixture{N_,F,S,T}}} ) where {N_,F<:AbstractRelativeRoots,S,T}
#   #
#   # return nothing

#   p = Int[1:ccwl.xDim;]
#   ccwl.partialDims = SVector(Int32.(p)...)

#   # should be done with constructor only 
#   for thrid in 1:Threads.nthreads()
#     # length(ccwl.cpt[thrid].p) != ccwl.xDim ? resize!(ccwl.cpt[thrid].p, ccwl.xDim) : nothing
#     ccwl.cpt[thrid].p = p  # SVector(Int32[1:ccwl.xDim;]...)
#   end
#   nothing
# end


##==============================================================================
## Deprecate code below before v0.24
##==============================================================================


@deprecate ContinuousMultivariate(N::Int) ContinuousEuclid(N)

# """
# $(TYPEDEF)

# Continuous variable of dimension `.dims` on manifold `.manifolds`.
# """
# struct ContinuousMultivariate{T1 <: Tuple} <: InferenceVariable
#   dims::Int
#   manifolds::T1
# end

# function ContinuousMultivariate(x::Int;
#                                 manifolds::T1=convert(Tuple, Euclidean(1))  )  where {T1 <: Tuple}
#   #
#   maniT = length(manifolds) < x ? ([manifolds[1] for i in 1:x]...,) : manifolds
#   ContinuousMultivariate{typeof(maniT)}(x, maniT)
# end

# # Legacy support
# getManifolds(::InstanceType{Manifolds.Euclidean{Tuple{N}, ℝ}}) where N = tuple([:Euclid for i in 1:N]...)
# getManifolds(::InstanceType{Manifolds.Circle{ℝ}})  = (:Circular,)


@deprecate solveFactorGraphParametric(w...; kw...) solveGraphParametric(w...; kw...)
@deprecate solveFactorGraphParametric!(fg::AbstractDFG; init::Bool=true, kwargs...) solveGraphParametric!(fg; init=init, kwargs...)



#
