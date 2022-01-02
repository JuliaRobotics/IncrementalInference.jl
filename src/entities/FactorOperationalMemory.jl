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
- Follow the Github project in IIF to better consolidate CCW FMD CPT CF CFM

Related 

[`CalcFactorMahalanobis`](@ref), [`CommonConvWrapper`](@ref), [`FactorMetadata`](@ref), [`ConvPerThread`](@ref)
"""
struct CalcFactor{T <: AbstractFactor, M, P <: Union{<:Tuple,Nothing,AbstractVector}, X}
  """ the interface compliant user object functor containing the data and logic """
  factor::T
  """ the metadata to be passed to the user residual function """
  metadata::M
  """ what is the sample (particle) id for which the residual is being calculated """
  _sampleIdx::Int
  """ legacy support when concerned with how many measurement tuple elements are used by user  """
  _measCount::Int
  """ legacy suport for measurement sample values of old functor residual functions """
  _legacyMeas::P
  """ legacy support for variable values old functor residual functions """
  _legacyParams::X
  """ allow threading for either sampling or residual calculations (workaround for thread yield issue) """
  _allowThreads::Bool
end


"""
    $TYPEDEF

Internal parametric extension to [`CalcFactor`](@ref) used for buffering measurement and calculating Mahalanobis distance

Related

[`CalcFactor`](@ref)
"""
struct CalcFactorMahalanobis{CF<:CalcFactor, S, N}
  calcfactor!::CF
  varOrder::Vector{Symbol}
  meas::NTuple{N, <:AbstractVector{Float64}}
  iΣ::NTuple{N, Matrix{Float64}}
  specialAlg::S
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
- Follow the Github project in IIF to better consolidate CCW FMD CPT CF CFM
- TODO standardize -- #927, #1025, #784, #692, #640
- TODO make immutable #825
"""
mutable struct FactorMetadata{FV<:AbstractVector{<:DFGVariable}, 
                              VL<:AbstractVector{Symbol}, 
                              AR<:NamedTuple, 
                              CD}
  # full list of Vector{DFGVariable} connected to the factor
  fullvariables::FV # Vector{<:DFGVariable}
  #TODO full variable can perhaps replace this
  variablelist::VL # Vector{Symbol} 
  # TODO rename/consolidate field in CCW,
  arrRef::AR
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
  # slight numerical perturbation for degenerate solver cases such as division by zero
  perturb::Vector{Float64}
  # working memory location for optimization routines on target decision variables
  X::Vector{P}
  # working memory to store residual for optimization routines
  res::R # was Vector{Float64}
end



"""
$(TYPEDEF)

Main factor memory container used during inference operations -- i.e. values specific to one complete convolution operation

Notes
- CCW does not get serialized / persisted
- At writing, the assumption is there is just one CCW per factor
- Any multithreaded design needs to happens as sub-constainers inside CCW or otherwise, to carry separate memory.
- Since #467, `CalcFactor` is the only type 'seen by user' during `getSample` or function residual calculations `(cf::CalcFactor{<:MyFactor})`, s.t. `MyFactor <: AbstractRelative___`
- There also exists a `CalcFactorMahalanobis` for parameteric computations using as much of the same mechanics as possible.

DevNotes
- Follow the Github project in IIF to better consolidate CCW FMD CPT CF CFM

Related 

[`CalcFactor`](@ref), [`CalcFactorMahalanobis`](@ref), [`FactorMetadata`](@ref), [`ConvPerThread`](@ref)
"""
mutable struct CommonConvWrapper{ T<:AbstractFactor,
                                  H<:Union{Nothing, Distributions.Categorical},
                                  C<:Union{Nothing, Vector{Int}},
                                  NTP <: NamedTuple,
                                  G,
                                  MT} <: FactorOperationalMemory
  #
  ### Values consistent across all threads during approx convolution
  usrfnc!::T # user factor / function
  # general setup
  xDim::Int
  zDim::Int
  # is this a partial constraint as defined by the existance of factor field `.partial::Tuple`
  partial::Bool 
  # multi hypothesis settings
  hypotheses::H 
  # categorical to select which hypothesis is being considered during convolution operation
  certainhypo::C
  nullhypo::Float64
  # parameters passed to each hypothesis evaluation event on user function, #1321
  params::NTP 
  # which index is being solved for in params?
  varidx::Int 
  # FIXME make type stable, JT should now be type stable if rest works
  #   user defined measurement values for each approxConv operation
  measurement::Vector{MT} 
  # TODO refactor and deprecate this old approach, Union{Type{SingleThreaded}, Type{MultiThreaded}}
  threadmodel::Type{<:_AbstractThreadModel} 
  # FIXME, deprecate for only `readonly` and build CalcFactor objects on stack instead
    ## will be obsolete: particular convolution computation values per particle idx (varies by thread)
  cpt::Vector{<:ConvPerThread}
  # inflationSpread
  inflation::Float64
  # Which dimensions does this factor influence.  Sensitive (mutable) to both which 'solvefor index' variable and whether the factor is partial dimension
  partialDims::Vector{<:Integer}
  # variable types for points in params
  vartypes::Vector{DataType}
  # experimental feature to embed gradient calcs with ccw
  _gradients::G
end



#
