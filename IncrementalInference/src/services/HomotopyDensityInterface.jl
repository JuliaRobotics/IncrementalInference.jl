

function setBelief!(
  state::State, 
  hode::ApproxManifoldProducts.HomotopyDensity, 
  setinit::Bool=true, 
  # ipc::AbstractVector{<:Real}=[0.0;]; # TODO, remove, use hode.observability, drop separate ipc here
  # solveKey::Symbol = :default
)
  @assert getStateKind(state) == getStateKind(hode) "statekind (i.e. lazy manifold serde) mismatch between variable and incoming belief $(getStateKind(vari)) vs $(getStateKind(hode))"
  state.belief = convert(typeof(state.belief), hode)
  setinit ? (state.initialized = true) : nothing
  nothing
end


function setBelief!(
  vari::VariableCompute, 
  hode::ApproxManifoldProducts.HomotopyDensity, 
  setinit::Bool=true, 
  ipc::AbstractVector{<:Real}=[0.0;]; # TODO, remove, use hode.observability, drop separate ipc here
  solveKey::Symbol = :default
)
  state = getState(vari, solveKey)
  setBelief!(state, hode, setinit)
end


## ==============================================================================================
## MOVE TO / CONSOLIDATE WITH DFG
## ==============================================================================================

# """
#     $(SIGNATURES)

# Fetch the variable marginal joint sampled points.  Use [`getBelief`](@ref) to retrieve the full Belief object.
# """
# #FIXME replace with refPoints
getVal(v::VariableCompute; solveKey::Symbol = :default) = getPoints(getState(v, solveKey).belief; permute=false)
function getVal(v::VariableCompute, idx::Int; solveKey::Symbol = :default)
  return getPoints(getState(v, solveKey).belief; permute=false)[idx]
end
getVal(vnd::State) = getPoints(vnd; permute=false)
getVal(vnd::State, idx::Int) = getPoints(vnd; permute=false)[idx]
function getVal(dfg::AbstractDFG, lbl::Symbol; solveKey::Symbol = :default)
  return getVal(getVariable(dfg, lbl); solveKey)
end

"""
    $(SIGNATURES)

Get the number of points used for the current marginal belief estimate represtation for a particular variable in the factor graph.
"""
function getNumPts(v::VariableCompute; solveKey::Symbol = :default)::Int
  return Npts(getState(v, solveKey).belief)
end

# function AMP.getBW(vnd::State)
#   return ApproxManifoldProducts.getBW(vnd.belief)
# end

# # setVal! assumes you will update values to database separate, this used for local graph mods only
# function getBWVal(v::VariableCompute; solveKey::Symbol = :default)
#   return getBW(getState(v, solveKey))
# end
# function setBW!(vd::State, bw::Array{Float64, 2}; solveKey::Symbol = :default)
#   DFG.refBandwidth(vd) .= bw # FIXME
#   return nothing
# end
# function setBW!(v::VariableCompute, bw::Array{Float64, 2}; solveKey::Symbol = :default)
#   setBW!(getState(v, solveKey), bw)
#   return nothing
# end

# function setVal!(
#   vd::State, 
#   val::AbstractVector{P}; 
#   observability::AbstractVector{<:Real} = [0.0;]
# ) where {P}
#   points = DFG.refPoints(vd)
#   resize!(points, length(val))
#   points .= val

#   observability = DFG.refObservability(vd)
#   resize!(observability, length(observability))
#   observability .= observability
#   return nothing
# end
# function setVal!(
#   v::VariableCompute,
#   val::AbstractVector{P};
#   solveKey::Symbol = :default,
#   observability::AbstractVector{<:Real} = [0.0;],
# ) where {P}
#   setVal!(getState(v, solveKey), val; observability)
#   return nothing
# end
# function setVal!(
#   vd::State,
#   val::AbstractVector{P},
#   bw::AbstractMatrix{Float64};
#   observability::AbstractVector{<:Real} = [0.0;],
# ) where {P}
#   setVal!(vd, val; observability)
#   setBW!(vd, bw)
#   return nothing
# end
# function setVal!(
#   v::VariableCompute,
#   val::AbstractVector{P},
#   bw::AbstractMatrix{Float64};
#   solveKey::Symbol = :default,
#   observability::AbstractVector{<:Real} = [0.0;],
# ) where {P}
#   setVal!(v, val; solveKey, observability)
#   setBW!(v, bw; solveKey)
#   return nothing
# end
# function setVal!(
#   vd::State,
#   val::AbstractVector{P},
#   bw::AbstractVector{Float64};
#   observability::AbstractVector{<:Real} = [0.0;],
# ) where {P}
#   setVal!(vd, val, reshape(bw, length(bw), 1); observability)
#   return nothing
# end
# function setVal!(
#   v::VariableCompute,
#   val::AbstractVector{P},
#   bw::AbstractVector{Float64};
#   solveKey::Symbol = :default,
#   observability::AbstractVector{<:Real} = [0.0;],
# ) where {P}
#   setVal!(getState(v, solveKey), val, bw; observability)
#   return nothing
# end
# function setVal!(
#   dfg::AbstractDFG,
#   sym::Symbol,
#   val::AbstractVector{P};
#   solveKey::Symbol = :default,
#   observability::AbstractVector{<:Real} = [0.0;],
# ) where {P}
#   return setVal!(getVariable(dfg, sym), val; solveKey, observability)
# end

"""
    $SIGNATURES

Set the point centers and bandwidth parameters of a variable node, also set `isInitialized=true` if `setinit::Bool=true` (as per default).

Notes
- `initialized` is used for initial solve of factor graph where variables are not yet initialized.
- `inferdim` is used to identify if the initialized was only partial.
"""
function setValKDE!(
  state::State,
  pts::AbstractVector{P},
  bws::Vector{Float64},
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;],
) where {P}
  #
  # @info "setValKDE!" string(bws)
  hode = HomotopyDensity_legacy(
    getStateKind(state),
    pts;
    bw = diagm(bws),
    Observability = ipc
  )
  setBelief!(state, hode, setinit)
  return nothing
end

function setValKDE!(
  state::State,
  val::AbstractVector{P},
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;],
) where {P}
  # recover variableType information
  varType = getStateKind(state)
  hode = HomotopyDensity_legacy(varType, val; observability=ipc)
  # p = AMP.manikde!(varType, val)
  setValKDE!(state, hode, setinit, ipc)
  return nothing
end

function setValKDE!(
  v::VariableCompute,
  val::AbstractVector{P},
  bws::AbstractMatrix, # obsolete -- from when bw diags were packed as columns in a matrix
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;];
  solveKey::Symbol = :default,
) where {P}
  # FIXME, use HomotopyDensity directly here instead
  setValKDE!(getState(v, solveKey), val, diag(bws), setinit, ipc)

  return nothing
end

function setValKDE!(
  v::VariableCompute,
  val::AbstractVector{P},
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;];
  solveKey::Symbol = :default,
) where {P}
  vnd = getState(v, solveKey)
  # recover variableType information
  setValKDE!(vnd, val, setinit, ipc)
  return nothing
end
function setValKDE!(
  v::VariableCompute,
  em::TreeBelief,
  setinit::Bool = true;
  # inferdim::Union{Float32, Float64, Int32, Int64}=0;
  solveKey::Symbol = :default,
)
  #
  setValKDE!(v, em.val, em.bw, setinit, em.infoPerCoord; solveKey = solveKey)
  return nothing
end
function setValKDE!(
  v::VariableCompute,
  mkd::ApproxManifoldProducts.HomotopyDensity,
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;];
  solveKey::Symbol = :default,
)
  #
  # @error("TESTING setValKDE! ", solveKey, string(listStates(v)))
  setValKDE!(getState(v, solveKey), mkd, setinit, Float64.(ipc))
  return nothing
end
function setValKDE!(
  dfg::AbstractDFG,
  sym::Symbol,
  mkd::ApproxManifoldProducts.HomotopyDensity,
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;];
  solveKey::Symbol = :default,
)
  #
  setValKDE!(getVariable(dfg, sym), mkd, setinit, ipc; solveKey = solveKey)
  return nothing
end

