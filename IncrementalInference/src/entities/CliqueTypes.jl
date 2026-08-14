# Clique types

# this is a developmental type, will be standardized after conclusion of #1010
# TODO resolve type instability
const MsgRelativeType = Vector{
  NamedTuple{(:variables, :likelihood), Tuple{Vector{Symbol}, <:DFG.AbstractRelativeObservation}},
}

const MsgPriorType = Dict{Symbol, MsgPrior{<:ManifoldKernelDensity}}

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

The concrete subtype says **what the density is**, which is what the receiver needs in order to know
what to do with it — the distinction is not cosmetic, and getting it wrong is silent (the solve still
converges, to the wrong covariance):

| subtype | is | receiver must |
|---|---|---|
| [`LinearizedLikelihood`](@ref) | `p(z_subtree │ S)` — evidence from **inside** the sender's subtree | pull back, then [`fuse`](@ref) |
| [`LinearizedBelief`](@ref) | `p(x │ z_all)` — the sender's full posterior | pull back, subtract its own upward message ([`calcCavityprecision`](@ref)), then use |

A cavity — `p(z_outside │ S)`, the complement of a likelihood — would be a third subtype, used
directly with no subtraction.  Nothing sends one yet.

Every subtype must be pulled back onto the receiver's tangent space first, since cliques re-linearize
independently ([`pullback`](@ref)).
"""
abstract type LinearizedContent end

#TODO rename, see #1954, options:
# TangentSeparatorPrior, TangentMarginal, GaugedPrior, LinearGaugeLikelihood
# LinearSeparatorPrior, MarginalLikelihood, GaugedMarginalPrior, GaugedGaussian

"""
    $TYPEDEF

**Upward** payload: the eliminated separator system, `p(z_subtree │ S)`, anchored at `bundle.point`.

A likelihood — it summarizes only the sender's subtree, so the receiver **fuses** it.

$(TYPEDFIELDS)
"""
struct LinearizedLikelihood <: LinearizedContent
  """the separator marginal being sent, with its base point"""
  bundle::DensityBundlePoint
  """did this clique **or anything in its subtree** re-linearize this pass?  The receiver must re-fuse
  even when its own base point held, because its `Λ` sums its children's."""
  relinearized::Bool
  """is this clique — or anything below it — still stepping by more than `tol`?

  `relinearized`: that asks "is my cached Jacobian stale" (`relinearizeTol`, a work decision), 
  `moving` asks "has the solve stopped" (`tol`, a convergence decision)."""
  moving::Bool
end

# both flags default to "assume the worst", so a caller that does not track them cannot accidentally
# claim a clique is settled
LinearizedLikelihood(bundle::DensityBundlePoint) = LinearizedLikelihood(bundle, true, true)
LinearizedLikelihood(bundle::DensityBundlePoint, relinearized::Bool) =
  LinearizedLikelihood(bundle, relinearized, true)

"""
    $TYPEDEF

**Downward** payload: the sender's joint posterior over its own variables, plus the solved tangent mean.

A belief — it already contains the receiver's own upward contribution, so the
receiver must subtract that ([`calcCavityprecision`](@ref)) before using the precision.  The *mean* half
needs no such correction, since back-substitution conditions on it rather than fusing it.

`Δ` stays a separate vector rather than folding into `η = ΛΔ`: a descendant's mean is largely
null-space content whenever the graph has gauge freedom, and `η` would discard exactly that.

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
Belief message for message passing on the tree.  This should be considered an incomplete joint probility.

Notes:
- belief -> Dictionary of [`TreeBelief`](@ref)
- variableOrder -> Ordered variable id list of the seperators in cliqueLikelihood
- cliqueLikelihood -> marginal distribution (<: `SamplableBelief`) over clique seperators.
- Older names include: productFactor, Fnew, MsgPrior, LikelihoodMessage

DevNotes:
- Used by both nonparametric and parametric.
- Objective for parametric case: `MvNormal(μ=[:x0;:x2;:l5], Σ=[+ * *; * + *; * * +])`.
- Part of the consolidation effort, see #459.
- Better conditioning for joint structure in the works using deconvolution, see #579, #635.
  - TODO confirm why <: Singleton.

$(TYPEDFIELDS)
"""
mutable struct LikelihoodMessage{T <: MessageType} <: AbstractPriorObservation
  sender::NamedTuple{(:id, :step), Tuple{Int, Int}}
  status::CliqStatus
  belief::Dict{Symbol, TreeBelief} # will eventually be deprecated
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

A clique's reusable residual/Jacobian machinery, plus the buffers they write into.  Call it through
[`linearize!`](@ref).

Lifetime is the tree's: rebuilding the tree clears it.

$(TYPEDFIELDS)
"""
struct CliqueLinearizer{MT, CF, JF}
  """joint manifold of the clique's variables, so callers need not rebuild it"""
  M::MT
  """stacked residual `r(p)` over the clique's own factors"""
  costF!::CF
  """Jacobian `J(p)` of that residual, in the tangent basis at `p`"""
  jacF!::JF
  """`m×n` Jacobian buffer, overwritten by every [`linearize!`](@ref)"""
  J::Matrix{Float64}
  """length-`m` residual buffer, overwritten by every [`linearize!`](@ref)"""
  r::Vector{Float64}
  """the clique's layout object"""
  layout::JointLayout
end

"""
    $SIGNATURES

Linearize the clique's own factors at `p`, returning `(J, r)`.

Return **views onto reused buffers**, overwritten by the next call.
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

What [`eliminateCliqueFrontals!`](@ref) produced for one clique, retained for the downward pass and for later sweeps.

Two axes name the fields:
 - **whose factors** (this clique alone, or the whole subtree rooted here) and
 - **which variables**

| | own factors only | own **+** every descendant's |
|---|---|---|
| all clique variables | `cliquelinearization` | `Λ_subtree` |
| separators only | — | `separatorlikelihood` |
| frontals │ separators | — | `conditional` |

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct CliqueElimination
  """frontal coordinate indices, in clique frontal order"""
  idx_F::UnitRange{Int}
  """separator coordinate indices, in **partition** order (see `separatorlikelihood.labels`)"""
  idx_S::UnitRange{Int}
  """this clique's own factors only, **before** child fusion.

  `cliquelinearization.point` is the **linearization point**: where the Jacobian was evaluated, not where
  the clique currently sits.  `_isLinearizationCurrent` measures the distance between the two against
  `relinearizeTol` to decide whether the Jacobian can be reused, which is the whole caching mechanism."""
  cliquelinearization::DensityBundlePoint
  """joint precision over `[frontals; separators]` for the whole subtree — `cliquelinearization`'s with
  every descendant's upward message fused in.  The difference is the fusion, not the variable set.
  Coordinate layout of `idx_F`/`idx_S`.

  Only precision needed because elimination consumes the joint information vector."""
  Λ_subtree::Matrix{Float64}
  """the subtree's likelihood over the **separators alone**, `p(z_subtree │ S)` — what goes on the wire
  upward as a [`LinearizedLikelihood`](@ref).

  Needed for the downward pass, derivable from `Λ_subtree` but kept as elimination already provides it.
  `Λ_subtree` is the joint over `F ∪ S`, which is what carries the frontal covariance; this is the `S`
  block of it.
  """
  separatorlikelihood::DensityBundlePoint
  """the conditional `p(F│S)`, completed on the down pass as `ΔF = v - W·ΔS`.

  `W` is recomputable from `Λ_subtree`, but `v` is not."""
  conditional::GaussianConditional

  """did the call that produced this **skip** the Jacobian and reuse a previous `cliquelinearization`?
  Describes that one call, so it is read only on the path that just made it."""
  reusedlinearization::Bool = false
  """norm of the frontal step this clique last took, `Inf` until it has solved once.  Movement happens
  on the DOWN pass, after the up message has gone, so it is reported upward on the NEXT sweep."""
  laststep::Float64 = Inf
  """did this clique re-linearize on the pass that wrote this — **or** any clique below it?
  Staleness propagates upward, so the ROOT's flag is the fixed-point test for the whole loop."""
  relinearized::Bool = true
  """is this clique — or anything below it — still stepping by more than `tol`?  The convergence
  reduction, carried here so the DOWN pass can read what the UP pass computed."""
  moving::Bool = true
  """The clique's state between sweeps: where it sits and what it believes (joint over all clique variables), 
  as a **base point plus a tangent offset** rather than a materialized point.

  After a downward pass this is what that pass solved — `Λ_subtree` with the cavity folded into its
  separator block, and the total offset from `cliquelinearization.point`.  `Δ` is a total, not an
  increment: back-substitution solves against the conditional built at that base, so successive sweeps
  that reuse a linearization overwrite it rather than accumulate.

  Anchored at the **pre-step** point because that is the tangent space children pull back from."""
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

  # Linearized solve option: the clique's factorization, reused by the downward pass.
  elimination::Union{Nothing, CliqueElimination}
  # the clique's coordinate layout — labels in clique order and their frontal/separator blocking.
  #  THE clique's copy: everything below points at this object, so `===` is the coverage check.
  #  Depends only on the clique's variables, rebuilt only when when the tree is.
  cliquelayout::Union{Nothing, JointLayout}
  # the clique's residual/Jacobian machinery, reused across sweeps.
  #  It survives a re-linearization, which `elimination` does not — see [`CliqueLinearizer`](@ref).
  linearizer::Union{Nothing, CliqueLinearizer}
  # """parametric relinearization sweeps this clique ran.  Only meaningful on a ROOT, which is where
  # termination is decided; lives on the tree rather than the state machine container so it survives a
  # pass and can be reported by `solveTreeParametric!`."""
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
