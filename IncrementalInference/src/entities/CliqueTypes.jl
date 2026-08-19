# Clique types

# this is a developmental type, will be standardized after conclusion of #1010
# TODO resolve type instability
const MsgRelativeType = Vector{
  NamedTuple{(:variables, :likelihood), Tuple{Vector{Symbol}, <:DFG.AbstractRelativeObservation}},
}

const MsgPriorType = Dict{Symbol, MsgPrior{<:ApproxManifoldProducts.HomotopyDensity}}

"""
    $TYPEDEF

Internal development types used during consolidation.  Stores relative and prior information making up a joint likelihood 
message passed upward on the Bayes tree.
"""
mutable struct _MsgJointLikelihood
  relatives::IIF.MsgRelativeType
  priors::IIF.MsgPriorType
end

"""
    $TYPEDEF

Payload of a [`LinearizedMessage`](@ref).

| subtype | is | receiver must |
|---|---|---|
| [`LinearizedLikelihood`](@ref) | `p(z_subtree │ S)` — evidence from the sender's subtree | pull back, then [`fuse`](@ref) |
| [`LinearizedBelief`](@ref) | `p(x │ z_all)` — the sender's full posterior | pull back, subtract its own upward message ([`calcCavityprecision`](@ref)), then use |
"""
abstract type LinearizedContent end

#TODO rename, see #1954, options:
# TangentSeparatorPrior, TangentMarginal, GaugedPrior, LinearGaugeLikelihood
# LinearSeparatorPrior, MarginalLikelihood, GaugedMarginalPrior, GaugedGaussian

"""
    $TYPEDEF

**Upward** payload: the eliminated separator system `p(z_subtree │ S)`, anchored at `bundle.point`.

$(TYPEDFIELDS)
"""
struct LinearizedLikelihood <: LinearizedContent
  """the separator marginal being sent, with its base point"""
  bundle::DensityBundlePoint
  """did this clique or anything below it re-linearize this pass? Receiver must re-fuse even if its own base point held."""
  relinearized::Bool
  """is this clique or anything below it still stepping by more than `tol`?"""
  moving::Bool
end

# both flags default to "assume the worst", so a caller that does not track them cannot accidentally
# claim a clique is settled
LinearizedLikelihood(bundle::DensityBundlePoint) = LinearizedLikelihood(bundle, true, true)
LinearizedLikelihood(bundle::DensityBundlePoint, relinearized::Bool) =
  LinearizedLikelihood(bundle, relinearized, true)

"""
    $TYPEDEF

**Downward** payload: the sender's joint posterior plus the solved tangent mean.
Receiver must subtract its own upward contribution via [`calcCavityprecision`](@ref) before using the precision.
`Δ` is kept separate from `η` to preserve null-space (gauge) content.

$(TYPEDFIELDS)
"""
struct LinearizedBelief <: LinearizedContent
  """the sender's joint precision, with its base point"""
  bundle::DensityBundlePoint
  """solved tangent mean over `bundle.labels`, relative to `bundle.point`"""
  Δ::Vector{Float64}
end

LinearizedBelief(bundle::DensityBundlePoint, Δ::AbstractVector) =
  LinearizedBelief(bundle, Vector{Float64}(Δ))

"""
  $(TYPEDEF)
Belief message for message passing on the tree.  Incomplete joint probability.

Notes
- Used by both nonparametric and parametric.
- See #459, #1010 for consolidation status.

$(TYPEDFIELDS)
"""
mutable struct LikelihoodMessage{T <: MessageType} <: AbstractPriorObservation
  sender::NamedTuple{(:id, :step), Tuple{Int, Int}}
  status::CliqStatus
  belief::Dict{Symbol, HomotopyDensity} # TODO, will eventually be deprecated, use joint likelihood instead, also see #1010
  variableOrder::Vector{Symbol}
  cliqueLikelihood::Union{Nothing, SamplableBelief}  # TODO drop the Union
  msgType::T
  hasPriors::Bool
  # this is different from belief[].inferdim, as the total available infer dims remaining during down msgs -- see #910
  childSolvDims::Dict{Int, Float64}
  # calc differential factors for joint in the child clique
  jointmsg::_MsgJointLikelihood
  # diffJoints::Vector{NamedTuple{(:variables, :likelihood), Tuple{Vector{Symbol},DFG.AbstractRelativeObservation}}}
  # Linearized message, `nothing` for the nonparametric.
  linearized::Union{Nothing, LinearizedContent}
end

"""
    $TYPEDEF

Cache messages being passed on the tree, one container per clique.

Notes
- See model 2 (?) on IIF #674
"""
mutable struct MessageBuffer
  # up receive message buffer (multiple children, multiple messages)
  upRx::Dict{Int, LikelihoodMessage}
  # down receive message buffer (one parent)
  downRx::Union{Nothing, LikelihoodMessage}
  # RESERVED up outgoing message buffer (one parent)
  upTx::Union{Nothing, LikelihoodMessage}
  # RESERVED down outgoing message buffer (multiple children but one message)
  downTx::Union{Nothing, LikelihoodMessage}
end
MessageBuffer() = MessageBuffer(Dict{Int, LikelihoodMessage}(), nothing, nothing, nothing)

##==============================================================================
## BayesTreeNodeData
##==============================================================================

"""
    $TYPEDEF

A clique's reusable residual/Jacobian machinery. Call via [`linearize!`](@ref). Lifetime is the tree's.

$(TYPEDFIELDS)
"""
struct CliqueLinearizer{MT, CF, JF}
  """joint manifold of the clique's variables"""
  M::MT
  """stacked residual `r(p)` over the clique's own factors"""
  costF!::CF
  """Jacobian `J(p)` of that residual in the tangent basis at `p`"""
  jacF!::JF
  """`m×n` Jacobian buffer, overwritten by every [`linearize!`](@ref)"""
  J::Matrix{Float64}
  """length-`m` residual buffer, overwritten by every [`linearize!`](@ref)"""
  r::Vector{Float64}
  """the clique's layout"""
  layout::JointLayout
end

"""
    $SIGNATURES

Linearize at `p`, returning `(J, r)` as views onto reused buffers (overwritten on the next call).
"""
function linearize!(linearizer::CliqueLinearizer, p)
  linearizer.jacF!(linearizer.M, linearizer.J, p)
  linearizer.costF!(linearizer.M, linearizer.r, p)
  return (linearizer.J, linearizer.r)
end

covers(linearizer::Union{Nothing, CliqueLinearizer}, layout::JointLayout) =
  isnothing(linearizer) || linearizer.layout === layout

"""
    $TYPEDEF

What [`eliminateCliqueFrontals!`](@ref) produced for one clique, retained for the downward pass and later sweeps.

Fields span two axes — *whose factors* (own vs. whole subtree) and *which variables*:

| | own factors only | own + every descendant's |
|---|---|---|
| all clique variables | `cliquelinearization` | `Λ_subtree` |
| separators only | — | `separatorlikelihood` |
| frontals │ separators | — | `conditional` |

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct CliqueElimination
  """frontal coordinate indices"""
  idx_F::UnitRange{Int}
  """separator coordinate indices, in partition order (see `separatorlikelihood.labels`)"""
  idx_S::UnitRange{Int}
  """own factors only, before child fusion, anchored at the linearization point.
  `_isLinearizationCurrent` checks distance from this point to decide whether to reuse the Jacobian."""
  cliquelinearization::DensityBundlePoint
  """joint precision over `[frontals; separators]` for the whole subtree (own + all children fused in)"""
  Λ_subtree::Matrix{Float64}
  """subtree likelihood over separators `p(z_subtree │ S)` — the upward wire message"""
  separatorlikelihood::DensityBundlePoint
  """conditional `p(F│S)`; completed on the down pass as `ΔF = v - W·ΔS`. `v` is not rederivable from `Λ_subtree`."""
  conditional::GaussianConditional

  """was the Jacobian reused (not re-evaluated) for this elimination?"""
  reusedlinearization::Bool = false
  """norm of the last frontal step; `Inf` until the clique has solved once"""
  laststep::Float64 = Inf
  """did this clique or anything below it re-linearize this pass? Root's flag is the whole-tree fixed-point test."""
  relinearized::Bool = true
  """is this clique or anything below it still stepping by more than `tol`?"""
  moving::Bool = true
  """clique state between sweeps: base point + tangent offset over all clique variables.
  Written by the downward pass; `Δ` is a total (not an increment) anchored at `cliquelinearization.point`."""
  posterior::LinearizedBelief
end

"""
$(TYPEDEF)

Data structure for each clique in the Bayes (Junction) tree.
"""
mutable struct BayesTreeNodeData
  status::CliqStatus
  frontalIDs::Vector{Symbol}
  separatorIDs::Vector{Symbol}
  inmsgIDs::Vector{Symbol} # Int
  potIDs::Vector{Symbol} # Int # this is likely redundant TODO -- remove
  potentials::Vector{Symbol}
  partialpotential::Vector{Bool}

  dwnPotentials::Vector{Symbol}
  dwnPartialPotential::Vector{Bool}

  cliqAssocMat::Array{Bool, 2}
  cliqMsgMat::Array{Bool, 2}
  directvarIDs::Vector{Symbol}
  directFrtlMsgIDs::Vector{Symbol}
  msgskipIDs::Vector{Symbol}
  itervarIDs::Vector{Symbol}
  directPriorMsgIDs::Vector{Symbol}
  debug::Any
  debugDwn::Any

  allmarginalized::Bool
  initialized::Symbol
  upsolved::Bool
  downsolved::Bool
  isCliqReused::Bool             # holdover

  # JT Local messages saved for cache and debugging, see IIF #675
  messages::MessageBuffer

  # Linearized solve: clique factorization, reused by the downward pass
  elimination::Union{Nothing, CliqueElimination}
  # THE clique's layout; everything below points at this object (`===` is the coverage check)
  cliquelayout::Union{Nothing, JointLayout}
  # clique's residual/Jacobian machinery, survives re-linearization (unlike `elimination`)
  linearizer::Union{Nothing, CliqueLinearizer}
  parIter::Int
end

## Packed types for serialization

mutable struct PackedBayesTreeNodeData
  frontalIDs::Vector{Symbol}
  separatorIDs::Vector{Symbol}
  inmsgIDs::Vector{Symbol} # Int
  potIDs::Vector{Symbol} # Int # this is likely redundant TODO -- remove
  potentials::Vector{Symbol}
  partialpotential::Vector{Bool}
  dwnPotentials::Vector{Symbol}
  dwnPartialPotential::Vector{Bool}
  cliqAssocMat::Array{Bool, 2}
  cliqMsgMat::Array{Bool, 2}
  directvarIDs::Vector{Symbol} # Int
  directFrtlMsgIDs::Vector{Symbol} # Int
  msgskipIDs::Vector{Symbol} # Int
  itervarIDs::Vector{Symbol} # Int
  directPriorMsgIDs::Vector{Symbol} # Int
end

## Full Clique Types

struct CliqueId{T}
  value::T
end

"""
    $(TYPEDEF)
Structure to store clique data
DEV NOTES: To replace TreeClique completely
    $(FIELDS)
"""
mutable struct TreeClique
  "Interger id unique within a tree with userId, robotId, sessionId"
  id::CliqueId{Int64} # not to be confused with the underlying index used by LightGraphs.jl, see issue #540
  "Data as `BayesTreeNodeData`"
  data::BayesTreeNodeData
  "Drawing attributes"
  attributes::Dict{String, Any}
  #solveInProgress #on a clique level a "solve in progress" might be very handy
end

#
