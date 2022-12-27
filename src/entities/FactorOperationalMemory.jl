
"""
$TYPEDEF

User factor interface method for computing the residual values of factors.

Notes
- Also see #467 on API consolidation

```julia
function (cf::CalcFactor{<:LinearRelative})(res::AbstractVector{<:Real}, z, xi, xj)
  cf.variablelist
  cf.cache
  # generic on-manifold residual function 
  return distance(z, distance(xj, xi))
end
```

DevNotes
- Follow the Github project in IIF to better consolidate CCW FMD CPT CF CFM

Related 

[`CalcFactorMahalanobis`](@ref), [`CommonConvWrapper`](@ref)
"""
struct CalcFactor{FT <: AbstractFactor, X, C, VT <: Tuple}
  """ the interface compliant user object functor containing the data and logic """
  factor::FT
  """ what is the sample (particle) id for which the residual is being calculated """
  _sampleIdx::Int
  """ legacy support for variable values old functor residual functions.
      TBD, this is still being used by DERelative factors. """
  _legacyParams::X
  """ allow threading for either sampling or residual calculations (workaround for thread yield issue) """
  _allowThreads::Bool
  """ user cache of arbitrary type, overload the [`preambleCache`](@ref) function. NOT YET THREADSAFE """
  cache::C

  ## TODO Consolidation WIP with FactorMetadata
  # full list of variables connected to the factor
  # TODO make sure this list is of the active hypo only
  fullvariables::VT # Vector{<:DFGVariable} # FIXME change to tuple for better type stability
  # which index is being solved for?
  solvefor::Int
end

# should probably deprecate the abstract type approach?
abstract type _AbstractThreadModel end

"""
$(TYPEDEF)
"""
struct SingleThreaded <: _AbstractThreadModel end
"""
$(TYPEDEF)
"""
struct MultiThreaded <: _AbstractThreadModel end

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

[`CalcFactor`](@ref), [`CalcFactorMahalanobis`](@ref)
"""
mutable struct CommonConvWrapper{
  T <: AbstractFactor, 
  NTP <: Tuple, 
  G, 
  MT, 
  CT,
  HP <: Union{Nothing, <:Distributions.Categorical{Float64, Vector{Float64}}},
  CH <: Union{Nothing, Vector{Int}},
  VT <: Tuple
} <: FactorOperationalMemory
  #
  ### Values consistent across all threads during approx convolution
  usrfnc!::T # user factor / function
  # general setup
  xDim::Int
  zDim::Int
  # is this a partial constraint as defined by the existance of factor field `.partial::Tuple`
  partial::Bool
  # multi hypothesis settings #NOTE no need for a parameter as type is known from `parseusermultihypo`
  hypotheses::HP
  # categorical to select which hypothesis is being considered during convolution operation
  certainhypo::CH
  nullhypo::Float64
  # parameters passed to each hypothesis evaluation event on user function, #1321
  params::NTP # TODO rename to varValsLink::NTP
  # which index is being solved for in params?
  varidx::Int
  # FIXME make type stable, JT should now be type stable if rest works
  #   user defined measurement values for each approxConv operation
  measurement::Vector{MT}
  # # TODO refactor and deprecate this old approach, Union{Type{SingleThreaded}, Type{MultiThreaded}}
  # threadmodel::Type{<:_AbstractThreadModel}
  # inflationSpread
  inflation::Float64
  # Which dimensions does this factor influence.  Sensitive (mutable) to both which 'solvefor index' variable and whether the factor is partial dimension
  partialDims::Vector{<:Integer}
  # # variable types for points in params
  # vartypes::Vector{<:DataType}
  # experimental feature to embed gradient calcs with ccw
  _gradients::G
  # type used for cache
  dummyCache::CT

  #Consolidation from FMD
  fullvariables::VT # Vector{<:DFGVariable}

  #Consolidation from CPT
  # the actual particle being solved at this moment
  particleidx::Int
  # subsection indices to select which params should be used for this hypothesis evaluation
  activehypo::Vector{Int}
  # working memory to store residual for optimization routines
  res::Vector{Float64}
end

function Base.getproperty(ccw::CommonConvWrapper, f::Symbol)
  if f == :threadmodel
    return SingleThreaded
  elseif f == :vartypes
    return typeof.(getVariableType.(ccw.fullvariables))
  else
    return getfield(ccw, f)
  end
end

#
