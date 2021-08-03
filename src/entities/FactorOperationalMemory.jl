import Base: convert
import Base: ==

export CalcFactor


"""
$TYPEDEF

User factor interface method for computing the residual values of factors.

Notes
- Also see #467 on API consolidation

```juila
function (cf::CalcFactor{<:LinearRelative})(res::AbstractVector{<:Real}, z, xi, xj)
  cf.metadata.variablelist
  cf.metadata.targetvariable
  cf.metadata.usercache
  generic on-manifold residual function 
  
  return distance(z, distance(xj, xi))
end
```

DevNotes
- See IIF Project to consolidate CCW, CF, FMD, CPT

Related 

[`CommonConvWrapper`](@ref), [`FactorMetadata`](@ref), [`ConvPerThread`](@ref)
"""
struct CalcFactor{T <: AbstractFactor, M, P <: Union{<:Tuple,Nothing}, X <: AbstractVector}
  # the interface compliant user object functor containing the data and logic
  factor::T
  # the metadata to be passed to the user residual function
  metadata::M
  # what is the sample (particle) id for which the residual is being calculated
  _sampleIdx::Int
  # legacy support when concerned with how many measurement tuple elements are used by user 
  _measCount::Int
  # legacy suport for measurement sample values of old functor residual functions
  _legacyMeas::P
  # legacy support for variable values old functor residual functions
  _legacyParams::X
end



abstract type _AbstractThreadModel end

"""
$(TYPEDEF)
"""
struct SingleThreaded <: _AbstractThreadModel
end
"""
$(TYPEDEF)
"""
struct MultiThreaded <: _AbstractThreadModel
end



"""
$(TYPEDEF)

Notes
- type-stable `usercache`, see https://github.com/JuliaRobotics/IncrementalInference.jl/issues/783#issuecomment-665080114 

DevNotes
- TODO why not just a NamedTuple? Perhaps part of #467
- TODO better consolidate with CCW, CPT, CalcFactor
- TODO standardize -- #927, #1025, #784, #692, #640
- TODO make immutable #825
"""
mutable struct FactorMetadata{FV<:AbstractVector{<:DFGVariable}, 
                              VL<:AbstractVector{Symbol}, 
                              AR<:AbstractVector{<:AbstractArray}, 
                              CD}
  # full list of Vector{DFGVariable} connected to the factor
  fullvariables::FV # Vector{<:DFGVariable}
  #TODO full variable can perhaps replace this
  variablelist::VL # Vector{Symbol} 
  # TODO consolidate, same as ARR used in CCW,
  arrRef::AR # Vector{Matrix{Float64}}
  # label of which variable is being solved for
  solvefor::Symbol       
  # for type specific user data, see (? #784)
  cachedata::CD
end

"""
$(TYPEDEF)

DevNotes
- FIXME consolidate with CalcFactor
- TODO consolidate with CCW, FMd, CalcFactor
- TODO consider renaming `.p` to `.decisionDims`
  - TODO `.decisionDims::DD where DD <: Union{<:AbstractVector{Int},Colon}` -- ensure type stability
- TODO figure out if we want static parameter THRID
- TODO make static params {XDIM, ZDIM, P}
- TODO make immutable
"""
mutable struct ConvPerThread{R,F<:FactorMetadata,P}
  thrid_::Int
  # the actual particle being solved at this moment
  particleidx::Int
  # additional/optional data passed to user function
  factormetadata::F
  # subsection indices to select which params should be used for this hypothesis evaluation
  activehypo::Vector{Int}
  # Select which decision variables to include in a particular optimization run
  p::Vector{Int}
  # slight numerical perturbation for degenerate solver cases such as division by zero
  perturb::Vector{Float64}
  # working memory location for optimization routines on target decision variables
  X::Vector{P}
  # working memory to store residual for optimization routines
  res::R # was Vector{Float64}
end


function ConvPerThread( X::AbstractVector{P},
                        zDim::Int,
                        factormetadata::FactorMetadata;
                        particleidx::Int=1,
                        activehypo= 1:length(params),
                        p::AbstractVector{<:Integer}=collect(1:1),
                        perturb=zeros(zDim),
                        res=zeros(zDim),
                        thrid_ = 0  ) where P
  #
  return ConvPerThread{typeof(res), typeof(factormetadata), Any}( thrid_,
                        particleidx,
                        factormetadata,
                        Int[activehypo;],
                        Int[p...;],
                        perturb,
                        X,
                        res )
end



"""
$(TYPEDEF)
"""
mutable struct CommonConvWrapper{ T<:FunctorInferenceType,
                                  H<:Union{Nothing, Distributions.Categorical},
                                  C<:Union{Nothing, Vector{Int}},
                                  P,
                                  G} <: FactorOperationalMemory
  #
  ### Values consistent across all threads during approx convolution
  usrfnc!::T # user factor / function
  # general setup
  xDim::Int
  zDim::Int
  # special case settings
  specialzDim::Bool # is there a special zDim requirement -- defined by user
  # is this a partial constraint -- defined by user, see new 
  partial::Bool 
  # multi hypothesis settings
  hypotheses::H # categorical to select which hypothesis is being considered during convolution operation
  certainhypo::C
  nullhypo::Float64
  # values specific to one complete convolution operation
  # FIXME ? JT - What if all points are not on the same manifold?  See #1321
  #          DF, just make it NamedTuple? -- some detail on pinning CCW down at construction only
  params::Vector{<:AbstractVector{P}} # parameters passed to each hypothesis evaluation event on user function
  varidx::Int # which index is being solved for in params?
  # FIXME make type stable
  measurement::Tuple # user defined measurement values for each approxConv operation
  threadmodel::Type{<:_AbstractThreadModel} # Union{Type{SingleThreaded}, Type{MultiThreaded}}
  ### particular convolution computation values per particle idx (varies by thread)
  cpt::Vector{<:ConvPerThread}
  # inflationSpread
  inflation::Float64
  # DONT USE THIS YET which dimensions does this factor influence
  partialDims::Vector{Int} # should become SVector{N, Int}
  # variable types for points in params
  vartypes::Vector{DataType}
  # experimental feature to embed gradient calcs with ccw
  _gradients::G
end


function CommonConvWrapper( fnc::T,
                            X::AbstractVector{P},
                            zDim::Int,
                            params::AbstractVector{<:AbstractVector{Q}},
                            factormetadata::FactorMetadata;
                            specialzDim::Bool=false,
                            partial::Bool=false,
                            hypotheses::H=nothing,
                            certainhypo=nothing,
                            activehypo= 1:length(params),
                            nullhypo::Real=0,
                            varidx::Int=1,
                            measurement::Tuple=(Vector{Vector{Float64}}(),),  # FIXME should not be a Matrix
                            particleidx::Int=1,
                            xDim::Int=size(X,1),
                            partialDims::AbstractVector{<:Integer}=collect(1:size(X,1)), # TODO make this SVector, and name partialDims
                            perturb=zeros(zDim),
                            res::AbstractVector{<:Real}=zeros(zDim),
                            threadmodel::Type{<:_AbstractThreadModel}=MultiThreaded,
                            inflation::Real=3.0,
                            vartypes=typeof.(getVariableType.(factormetadata.fullvariables)),
                            gradients=nothing) where {T<:FunctorInferenceType,P,H,Q}
  #
  return  CommonConvWrapper(fnc,
                            xDim,
                            zDim,
                            specialzDim,
                            partial,
                            hypotheses,
                            certainhypo,
                            Float64(nullhypo),
                            params,
                            varidx,
                            measurement,
                            threadmodel,
                            (i->ConvPerThread(X, zDim,factormetadata, particleidx=particleidx,
                                              activehypo=activehypo, p=partialDims, 
                                              perturb=perturb, res=res )).(1:Threads.nthreads()),
                            inflation,
                            partialDims,  # SVector(Int32.()...)
                            DataType[vartypes...],
                            gradients)
end



#
