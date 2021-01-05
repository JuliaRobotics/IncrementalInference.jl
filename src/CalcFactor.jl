# New factor interface, something perhaps like this

export CalcFactor

# Also see #467 on API consolidation
# function (cf::CalcFactor{<:LinearRelative})(res::AbstractVector{<:Real}, z, xi, xj)
#   # cf.metadata.variablelist...
#   # cf.metadata.targetvariable
#   # cf.metadata.usercache
#   # generic on-manifold residual function 
#   res .= distance(z, distance(xj, xi))
# end

"""
    $TYPEDEF

New generation user factor interface method for computing the residual values of factors.

Notes
- Under development and still experimental.  Expected to become default method in IIF v0.20.0
"""
struct CalcFactor{T <: FunctorInferenceType, M, P <: Tuple, X <: AbstractVector}
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



"""
(cf::CalcFactor{T,M})( res::AbstractVector{<:Real}, meas..., params...)

Standard interface for calling a factor calculation, as if `cf.factor(residual, noise_process, parameters)`,
where factors are either library standard or user out-of-library factor definitions.  See documentation for
more details and tutorials on using your own factors (designed to be as easy as possible).

Notes
- These residual calculations use used to find non-Gaussian / multimodal (incl. PPE) and conventional Gaussian estimates. 
- `cf.legacyMeas == (measparams[1:cf._measCount]...,)`

Example
```julia
# TBD
```
"""
function (cf::CalcFactor)( res, measparams... ) #where {T<:FunctorInferenceType,M,P<:Tuple,X<:AbstractVector} {T,M,P,X}
  #
  # NOTE this is a legacy interface
  cf.factor(res, cf.metadata, cf._sampleIdx, cf._legacyMeas, cf._legacyParams...)
end



"""
    $SIGNATURES

Sample the factor stochastic model `N::Int` times and store the samples in the preallocated `ccw.measurement` container.

DevNotes
- Use in place operations where possible and remember `measurement` is a `::Tuple`.
- TODO only works on `.threadid()==1` at present, see #1094
"""
function freshSamples(cf::CalcFactor{<:FunctorInferenceType}, 
                      N::Int=1  )
  #
  getSample(cf, N)
end





function Base.show(io::IO, x::CalcFactor)
  println(io, )
  printstyled(io, " CalcFactor:\n", color=:blue)
  println(io, "  .factor: ", typeof(x.factor))
end

Base.show(io::IO, ::MIME"text/plain", x::CalcFactor) = show(io, x)
#