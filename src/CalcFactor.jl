# New factor interface, something perhaps like this

export CalcFactor

# Also see #467 on API consolidation
# function (cf::CalcFactor{<:LinearRelative})(res::Vector, z, xi, xj)
#   # cf.metadata.variablelist...
#   # cf.metadata.targetvariable
#   # cf.metadata.usercache
#   # generic on-manifold residual function 
#   res .= distance(z, distance(xj, xi))
# end


struct CalcFactor{T <: FunctorInferenceType, M, P <: Tuple, X <: AbstractVector}
  factor::T
  metadata::M
  _sampleIdx::Int
  _measCount::Int
  _legacyMeas::P
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
function (cf::CalcFactor{T,M})( res::AbstractVector{<:Real},
                                measparams... ) where {T<:FunctorInferenceType,M}
  #
  # NOTE this is a legacy interface
  cf.factor(res, cf.metadata, cf._sampleIdx, cf._legacyMeas, cf._legacyParams...)
end






#