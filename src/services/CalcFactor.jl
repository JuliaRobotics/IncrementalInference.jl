# New factor interface, something perhaps like this


export calcFactorResidualTemporary


# Helper function to construct CF from a CCW
CalcFactor(ccwl::CommonConvWrapper) = CalcFactor( ccwl.usrfnc!, _getFMdThread(ccwl), 0, length(ccwl.measurement), ccwl.measurement, ccwl.params)



"""
    $SIGNATURES

Sample the factor stochastic model `N::Int` times and store the samples in the preallocated `ccw.measurement` container.

DevNotes
- Use in place operations where possible and remember `measurement` is a `::Tuple`.
- TODO only works on `.threadid()==1` at present, see #1094
- Also see, JuliaRobotics/RoME.jl#465
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
residual = calcFactorResidualTemporary(Pose2Pose2(...), (RoME.Pose2,RoME.Pose2), (z_i,), (x1, x2))
```

Related

[`calcFactorResidual`](@ref), [`CalcResidual`](@ref), [`_evalFactorTemporary!`](@ref), [`approxConv`](@ref), [`_buildGraphByFactorAndTypes!`](@ref)
"""
function calcFactorResidualTemporary( fct::AbstractRelative, 
                                      varTypes::Tuple,
                                      measurement::Tuple,
                                      pts::Tuple;
                                      tfg::AbstractDFG = initfg(),
                                      _blockRecursion::Bool=false )
  #

  # build a new temporary graph
  _, _dfgfct = _buildGraphByFactorAndTypes!(fct, varTypes, pts, dfg=tfg, _blockRecursion=_blockRecursion)
  
  # get a fresh measurement if needed
  _measurement = if length(measurement) != 0
    measurement
  else
    # now use the CommonConvWrapper object in `_dfgfct`
    ccw = IIF._getCCW(_dfgfct)
    cfo = CalcFactor(ccw)
    getSample(cfo, 1)
  end

  # assume a single sample point is being run
  return calcFactorResidual(_dfgfct, _measurement..., pts...)
end


## =============================================================================================
## FactorOperationalMemory helper constructors
## =============================================================================================


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