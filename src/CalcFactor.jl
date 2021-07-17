# New factor interface, something perhaps like this

export CalcFactor
export testFactorResidualBinary

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

Evaluate the residual function for a single sample.

Notes
- Binary factors only at this stage, and `multihypo` does not have to be considered in this test
- Assumes calculation is for a single particle, so `meas::Tuple{Z,other}` is only a single particles value.

DevNotes
- TODO generalize for n-ary factors

Example
```julia
residual = testFactorResidualBinary(Pose2Pose2(...), (z_i,), (RoME.Pose2, x1), (RoME.Pose2, x2))
```

Related

[`approxConv`](@ref), [`CalcResidual`](@ref), [`_evalFactorTemporary!`](@ref)
"""
function testFactorResidualBinary(fct::AbstractRelative, 
                                  meas::Tuple,
                                  T_param_args... )
  #

  # TODO generalize beyond binary
  T1 = T_param_args[1][1]
  param1 = T_param_args[1][2]
  T2 = T_param_args[2][1]
  param2 = T_param_args[2][2]

  fg_ = initfg()
  X0 = addVariable!(fg_, :x0, T1)
  X1 = addVariable!(fg_, :x1, T2)
  addFactor!(fg_, [:x0;:x1], fct, graphinit=false)

  ccw = IIF._getCCW(fg_, :x0x1f1)
  
  ARR = getPoints.(getBelief.(fg_, [:x0;:x1]))
  fmd = FactorMetadata([X0;X1],[:x0;:x1],ARR,:x1,nothing)
  cfo = CalcFactor(fct, fmd, 1, 1, meas, ARR)

  # get a fresh measurement if needed
  meas = length(meas) != 0 ? meas : getSample(cfo, 1)
  cfo = CalcFactor(fct, fmd, 1, 1, meas, ARR)

  # residual vector
  zdim = IIF._getZDim(ccw)
  res = zeros(zdim)

  # calc the residual
  @time res = cfo(meas..., param1, param2)

  return res
end




#