
"""
    CliqStatus
Clique status message enumerated type with status:
NULL, INITIALIZED, UPSOLVED, MARGINALIZED, DOWNSOLVED, UPRECYCLED, ERROR_STATUS
"""
@enum CliqStatus NULL INITIALIZED UPSOLVED MARGINALIZED DOWNSOLVED UPRECYCLED ERROR_STATUS


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
  softtype::T
  # TODO -- DEPRECATE
  manifolds::Tuple{Vararg{Symbol}}# TODO #459
end
TreeBelief(p::BallTreeDensity,
           inferdim::Real=0.0,
           softtype::T=ContinuousScalar(),
           manifolds=getManifolds(softtype)) where {T <: InferenceVariable} = TreeBelief{T}(getPoints(p), getBW(p), inferdim, softtype, manifolds)

TreeBelief(val::Array{Float64,2},
           bw::Array{Float64,2},
           inferdim::Real=0.0,
           softtype::T=ContinuousScalar(),
           manifolds=getManifolds(softtype)) where {T <: InferenceVariable} = TreeBelief{T}(val, bw, inferdim, softtype, manifolds)

function TreeBelief(vnd::VariableNodeData)
  TreeBelief( vnd.val, vnd.bw, vnd.inferdim, getSofttype(vnd), getManifolds(vnd) )
end

TreeBelief(vari::DFGVariable, solveKey=:default) = TreeBelief(getSolverData(vari, solveKey))

getManifolds(treeb::TreeBelief) = getManifolds(treeb.softtype)


"""
  $(TYPEDEF)
Belief message for message passing on the tree.

Notes:
- belief -> Dictionary of [`TreeBelief`](@ref)
- variableOrder -> Ordered variable id list of the seperators in cliqueLikelihood
- cliqueLikelihood -> marginal distribution (<: `SamplableBelief`) over clique seperators.

DevNotes:
- Objective for parametric: `MvNormal(μ=[:x0;:x2;:l5], Σ=[+ * *; * + *; * * +])`
- TODO confirm why <: Singleton
- #459
  $(TYPEDFIELDS)
"""
mutable struct LikelihoodMessage <: Singleton
  status::CliqStatus
  belief::Dict{Symbol, TreeBelief}
  variableOrder::Vector{Symbol}
  cliqueLikelihood::Union{Nothing,SamplableBelief}
end

# EARLIER NAMES INCLUDE: productFactor, Fnew, MsgPrior, LikelihoodMessage

LikelihoodMessage(status::CliqStatus) =
        LikelihoodMessage(status, Dict{Symbol, TreeBelief}(), Symbol[], nothing)

LikelihoodMessage(status::CliqStatus, varOrder::Vector{Symbol}, cliqueLikelihood::SamplableBelief) =
        LikelihoodMessage(status, Dict{Symbol, TreeBelief}(), varOrder, cliqueLikelihood)

LikelihoodMessage(;status::CliqStatus=NULL,
                   beliefDict::Dict=Dict{Symbol, TreeBelief}(),
                   variableOrder=Symbol[],
                   cliqueLikelihood=nothing ) =
        LikelihoodMessage(status, beliefDict, variableOrder, cliqueLikelihood)
#


# FIXME, better standardize intermediate types
# used during nonparametric CK preparation, when information from multiple siblings must be shared together
const IntermediateSiblingMessages = Vector{Tuple{BallTreeDensity,Float64}}
const IntermediateMultiSiblingMessages = Dict{Symbol, IntermediateSiblingMessages}


### EVERYTHING BELOW IS/SHOULD BE DEPRECATED


# TODO this is casing problems between nonparametric and parametric
# const BeliefMessage = LikelihoodMessage


# Deprecated, replaced by LikelihoodMessage
# TODO - remove

# Dict{Symbol,   -- is for variable label
#  Vector{       -- multiple msgs for the same variable
#   Symbol,      -- Clique index
#   Int,         -- Depth in tree
#   BTD          -- Belief estimate
#   inferredDim  -- Information count
#  }
const TempUpMsgPlotting = Dict{Symbol,Vector{Tuple{Symbol, Int, BallTreeDensity, Float64}}}





#
