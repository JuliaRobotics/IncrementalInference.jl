# New factor interface, something perhaps like this

export CalcFactor
export calcFactorResidualTemporary

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

User factor interface method for computing the residual values of factors.
"""
struct CalcFactor{T <: AbstractFactor, M, P <: Tuple, X <: AbstractVector}
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


CalcFactor(ccwl::CommonConvWrapper) = CalcFactor( ccwl.usrfnc!, _getFMdThread(ccwl), 0, length(ccwl.measurement), ccwl.measurement, ccwl.params)



"""
    $SIGNATURES

Sample the factor stochastic model `N::Int` times and store the samples in the preallocated `ccw.measurement` container.

DevNotes
- Use in place operations where possible and remember `measurement` is a `::Tuple`.
- TODO only works on `.threadid()==1` at present, see #1094
"""
function sampleFactor(cf::CalcFactor{<:AbstractFactor}, 
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



"""
    $SIGNATURES

Helper function for evaluating factor residual functions, by adding necessary `CalcFactor` wrapper.
  
Notes
- Factor must already be in a factor graph to work
- Will not yet properly support all multihypo nuances, more a function for testing
- Useful for debugging a factor. 

Example
```julia
fg = generateCanonicalFG_Kaess()

residual = calcFactorResidual(fg, :x1x2f1, [1.0], [0.0], [0.0])
```

Related

[`calcFactorResidualTemporary`](@ref), [`_evalFactorTemporary!`](@ref), [`evalFactor`](@ref), [`approxConv`](@ref)
"""
calcFactorResidual(dfgfct::DFGFactor, args...; ccw::CommonConvWrapper=IIF._getCCW(dfgfct)) = CalcFactor(ccw)(args...)
calcFactorResidual(dfg::AbstractDFG, fctsym::Symbol, args...) = calcFactorResidual(getFactor(dfg, fctsym), args...)


"""
    $SIGNATURES

Evaluate the residual function for a single sample.

Notes
- Binary factors only at this stage, and `multihypo` does not have to be considered in this test
- Assumes calculation is for a single particle, so `meas::Tuple{Z,other}` is only a single particles value.

DevNotes
- TODO generalize for n-ary factors

Example
```julia
residual = calcFactorResidualTemporary(Pose2Pose2(...), (z_i,), (RoME.Pose2, x1), (RoME.Pose2, x2))
```

Related

[`calcFactorResidual`](@ref), [`CalcResidual`](@ref), [`_evalFactorTemporary!`](@ref), [`approxConv`](@ref), [`_buildGraphByFactorAndTypes!`](@ref)
"""
function calcFactorResidualTemporary( fct::AbstractRelative, 
                                      measurement::Tuple,
                                      varTypes::Tuple,
                                      pts::Tuple;
                                      tfg::AbstractDFG = initfg() )
  #

  # build a new temporary graph
  _, _dfgfct = _buildGraphByFactorAndTypes!(fct, varTypes, pts, dfg=tfg)
  
  # get a fresh measurement if needed
  measurement = if length(measurement) != 0
    measurement
  else
    # now use the CommonConvWrapper object in `_dfgfct`
    ccw = IIF._getCCW(_dfgfct)
    cfo = CalcFactor(ccw)
    getSample(cfo, 1)
  end

  # assume a single sample point is being run
  return calcFactorResidual(_dfgfct, measurement..., pts...)
end



#