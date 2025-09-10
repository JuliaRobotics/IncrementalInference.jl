
import DistributedFactorGraphs: getVariableType

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

const SamplableBelief = Union{
  <:Distributions.Distribution,
  <:KDE.BallTreeDensity, # FIXME deprecate
  <:AMP.ManifoldKernelDensity,
  <:AliasingScalarSampler,
  <:FluxModelsDistribution,
  <:HeatmapGridDensity,
  <:LevelSetGridNormal,
}

#Supported types for parametric
const ParametricTypes = Union{Normal, MvNormal}

"""
    $TYPEDEF

INTERMEDIATE DATA STRUCTURE DURING REFACTORING.

Representation of the belief of a single variable.

Notes:
- we want to send the joint, this is just to resolve consolidation #459 first.
- Long term objective is single joint definition, likely called `LikelihoodMessage`.
"""
struct TreeBelief{T <: StateType, P, M <: MB.AbstractManifold}
  val::Vector{P}
  bw::Array{Float64, 2}
  infoPerCoord::Vector{Float64}
  # see DFG #603, variableType defines the domain and manifold as well as group operations for a variable in the factor graph
  variableType::T
  # TODO -- DEPRECATE
  manifold::M # Tuple{Vararg{Symbol}} # NOTE added during #459 effort
  # only populated during up as solvableDims for each variable in clique, #910
  solvableDim::Float64
end

function TreeBelief(
  p::ManifoldKernelDensity,
  ipc::AbstractVector{<:Real} = [0.0;],
  variableType::T = ContinuousScalar(),
  manifold = getManifold(variableType),
  solvableDim::Real = 0,
) where {T <: StateType}
  return TreeBelief(getPoints(p), getBW(p), ipc, variableType, manifold, solvableDim)
end

function TreeBelief(
  val::AbstractVector{P},
  bw::Array{Float64, 2},
  ipc::AbstractVector{<:Real} = [0.0;],
  variableType::T = ContinuousScalar(),
  manifold::M = getManifold(variableType),
  solvableDim::Real = 0,
) where {P, T <: StateType, M <: MB.AbstractManifold}
  return TreeBelief{T, P, M}(val, bw, ipc, variableType, manifold, solvableDim)
end

function TreeBelief(vnd::State, solvDim::Real = 0)
  return TreeBelief(
    vnd.val,
    vnd.bw,
    vnd.infoPerCoord,
    getVariableType(vnd),
    getManifold(vnd),
    solvDim,
  )
end

function TreeBelief(vari::VariableCompute, solveKey::Symbol = :default; solvableDim::Real = 0)
  return TreeBelief(getState(vari, solveKey), solvableDim)
end
#

getVariableType(tb::TreeBelief) = tb.variableType

getManifold(treeb::TreeBelief) = getManifold(treeb.variableType)

function compare(t1::TreeBelief, t2::TreeBelief)
  TP = true
  TP = TP && norm(t1.val - t2.val) < 1e-5
  TP = TP && norm(t1.bw - t2.bw) < 1e-5
  TP = TP && isapprox(t1.infoPerCoord, t2.infoPerCoord; atol = 1e-4)
  TP = TP && t1.variableType == t2.variableType
  TP = TP && abs(t1.solvableDim - t2.solvableDim) < 1e-5
  return TP
end

#
