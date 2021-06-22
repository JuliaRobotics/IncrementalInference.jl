
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
## Deprecate code below before v0.27
##==============================================================================


# function manikde!(ptsArr::Vector{Vector{Float64}}, manis::Tuple)
#   arr = Matrix{Float64}(undef, length(ptsArr[1]), length(ptsArr))
#   @cast arr[i,j] = ptsArr[j][i]
#   manikde!(arr, manis)
# end


##==============================================================================
## Deprecate code below before v0.26
##==============================================================================

function getVal(vA::Vector{<:DFGVariable}, solveKey::Symbol=:default)
  @error "getVal(::Vector{DFGVariable}) is obsolete, use getVal.(DFGVariable) instead."
  # len = length(vA)
  # vals = Array{Array{Float64,2},1}()
  # cols = Array{Int,1}()
  # push!(cols,0)
  # rows = Array{Int,1}()
  # for v in vA
  #     push!(vals, getVal(v, solveKey=solveKey))
  #     c = size(vals[end],2)
  #     r = size(vals[end],1)
  #     push!(cols, floor(Int,c))
  #     push!(rows, floor(Int,r))
  # end
  # cols = cumsum(cols)
  # sc = cols[end]
  # rw = floor(Int,rows[1])
  # val = Array{Float64,2}(undef,rw, sc)
  # for i in 1:(len-1)
  #     val[:,(cols[i]+1):cols[i+1]] = vals[i]
  # end
  # val[:,(cols[len]+1):cols[len+1]] = vals[len] # and the last one
  # return val
end

# """
#     $SIGNATURES

# Get graph node (variable or factor) dimension.
# """
# DFG.getDimension(vartype::InferenceVariable) = vartype.dims #TODO Deprecate
# DFG.getDimension(vartype::Type{<:InferenceVariable}) = getDimension(vartype())


@deprecate ensureAllInitialized!(w...;kw...) initAll!(w...;kw...)

@deprecate getFactorMean(w...) IIF.getParametricMeasurement(w...)[1]

# """
#     $SIGNATURES

# Recover the mean (Gaussian) or estimate stochastic mean (non-Gaussian) value stored in a factor measurement.

# Related

# accumulateFactorMeans, solveBinaryFactorParameteric
# """
# function getFactorMean(fct::FunctorInferenceType)
#   fctt = typeof(fct)
#   error("no getFactorMean defined for $(fctt.name), has fields $(fieldnames(fctt))")
# end

# getFactorMean(fct::Normal) = [fct.μ]
# getFactorMean(fct::MvNormal) = fct.μ
# getFactorMean(fct::Union{<:BallTreeDensity,<:ManifoldKernelDensity}) = getKDEMean(fct)
# getFactorMean(fct::AliasingScalarSampler) = Statistics.mean(rand(fct,1000))

# getFactorMean(fct::DFGFactor) = getFactorMean(getFactorType(fct))

# getFactorMean(dfg::AbstractDFG, fctsym::Symbol) = getFactorMean(getFactor(dfg, fctsym))


# AMP.getManifolds(::T) where {T <: InferenceVariable} = getManifolds(getManifold(T))
# AMP.getManifolds(::Type{T}) where {T <: InferenceVariable} = getManifolds(getManifold(T))


# getManifolds(getManifold(T))
getManifolds(::InstanceType{T}) where {T <: Union{InferenceVariable, AbstractFactor}} = error("getManifolds is obsolete, use getManifold(...)::MB.AbstactManifold instead.")





#
