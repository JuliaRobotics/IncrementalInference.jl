
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


const SamplableBelief = Union{Distributions.Distribution, KDE.BallTreeDensity, AMP.ManifoldKernelDensity, AliasingScalarSampler, FluxModelsDistribution, HeatmapDensityRegular}

abstract type PackedSamplableBelief end

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
struct TreeBelief{T <: InferenceVariable}
  val::Array{Float64,2}
  bw::Array{Float64,2}
  inferdim::Float64
  # see DFG #603, variableType defines the domain and manifold as well as group operations for a variable in the factor graph
  variableType::T
  # TODO -- DEPRECATE
  manifolds::Tuple{Vararg{Symbol}} # NOTE added during #459 effort
  # only populated during up as solvableDims for each variable in clique, #910
  solvableDim::Float64 
end
TreeBelief( p::Union{<:BallTreeDensity, <:ManifoldKernelDensity},
            inferdim::Real=0,
            variableType::T=ContinuousScalar(),
            manifolds=getManifolds(variableType),
            solvableDim::Real=0) where {T <: InferenceVariable} = TreeBelief{T}(getPoints(p), getBW(p), inferdim, variableType, manifolds, solvableDim)

TreeBelief( val::Array{Float64,2},
            bw::Array{Float64,2},
            inferdim::Real=0,
            variableType::T=ContinuousScalar(),
            manifolds=getManifolds(variableType),
            solvableDim::Real=0) where {T <: InferenceVariable} = TreeBelief{T}(val, bw, inferdim, variableType, manifolds, solvableDim)

function TreeBelief(vnd::VariableNodeData, solvDim::Real=0)
  TreeBelief( vnd.val, vnd.bw, vnd.inferdim, getVariableType(vnd), getManifolds(vnd), solvDim )
end

TreeBelief(vari::DFGVariable, solveKey::Symbol=:default; solvableDim::Real=0) = TreeBelief( getSolverData(vari, solveKey) , solvableDim)

getVariableType(tb::TreeBelief) = tb.variableType

getManifolds(treeb::TreeBelief) = getManifolds(treeb.variableType)

function compare(t1::TreeBelief, t2::TreeBelief)
  TP = true
  TP = TP && norm(t1.val - t2.val) < 1e-5
  TP = TP && norm(t1.bw - t2.bw) < 1e-5
  TP = TP && abs(t1.inferdim - t2.inferdim) < 1e-5
  TP = TP && t1.variableType == t2.variableType
  TP = TP && abs(t1.solvableDim - t2.solvableDim) < 1e-5
  return TP
end





#
