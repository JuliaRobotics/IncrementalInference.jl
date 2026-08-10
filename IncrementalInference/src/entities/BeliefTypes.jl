


"""
    CliqStatus
Clique status message enumerated type with status.
"""
@enum CliqStatus NULL NO_INIT INITIALIZED UPSOLVED MARGINALIZED DOWNSOLVED UPRECYCLED ERROR_STATUS

# Used for UPWARD_DIFFERENTIAL, UPWARD_COMMON, DOWNWARD_COMMON marginalized types
abstract type MessagePassDirection end
struct UpwardPass <: MessagePassDirection end
struct DownwardPass <: MessagePassDirection end

abstract type MessageType end
struct NonparametricMessage <: MessageType end
struct ParametricMessage <: MessageType end

using DistributedFactorGraphs: PackedBelief

#TODO deprecate SamplableBelief
const SamplableBelief = Union{
  <:Distributions.Distribution,
  <:ApproxManifoldProducts.HomotopyDensity,
  <:AliasingScalarSampler,
  <:FluxModelsDistribution,
  <:HeatmapGridDensity,
  <:LevelSetGridNormal,
  <:Mixture,
}

#Supported types for parametric
const ParametricTypes = Union{Normal, MvNormal}

"""
    $TYPEDEF

INTERMEDIATE DATA STRUCTURE DURING REFACTORING.

Representation of the belief of a single variable.

DevNotes:
- we want to send the joint, this is just to resolve consolidation #459 first.
- Long term objective is single joint definition, likely called `LikelihoodMessage`.
- See 1929; wholesale replacement w `HomotopyDensity` over all clique dimensions.
"""
function TreeBelief(
  val::AbstractVector{P},
  bw::AbstractMatrix{Float64},
  ipc::AbstractVector{<:Real} = [0.0;],
  variableType::T = ContinuousScalar(),
  manifold::M = getManifold(variableType),
  solvableDim::Real = 0,
) where {P, T <: StateType, M <: MB.AbstractManifold}
  return HomotopyDensity_legacy(variableType, val; bw, observability = ipc)
end


function TreeBelief(
  p::ApproxManifoldProducts.HomotopyDensity,
  ipc::AbstractVector{<:Real} = [0.0;],
  variableType::T = ContinuousScalar(),
  manifold = getManifold(variableType),
  solvableDim::Real = 0,
) where {T <: StateType}
  return TreeBelief(getPoints(p), getBW(p), ipc, variableType, manifold, solvableDim)
end


function TreeBelief(state::State, solvDim::Real = 0)
  pts = getPoints(state.belief; permute=false) # TODO likely want to go back to sorted order here, DX debugging with permute=false
  cv = if length(getBW(state.belief)) <= 0 || !isassigned(getBW(state.belief), 1)
    @warn "TreeBelief constructor: belief has no bandwidth, defaulting to identity, incorrect -- must refactor, see #1929" maxlog = 20
    cov(state.belief)
  else
    getBW(state.belief)[1] # FIXME, bw for nonparametric, cov for parametric -- can reuse bw for parametric during AMP 15 refactor
  end
  obsv = DFG.refObservability(state)
  statekind = getStateKind(state)

  # @info "TreeBelief" string(pts[1]) string(cv) string(obsv)

  return TreeBelief(
    pts,
    cv,
    obsv,
    statekind,
    # getManifold(statekind),
    # solvDim,
  )
  # TreeBelief(DFG.getTopologyKind(state), state, solvDim)
end

function TreeBelief(vari::VariableCompute, solveKey::Symbol = :default; solvableDim::Real = 0)
  return TreeBelief(getState(vari, solveKey), solvableDim)
end
#



#
