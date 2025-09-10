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
