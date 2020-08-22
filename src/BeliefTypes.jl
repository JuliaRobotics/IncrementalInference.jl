
# """
#     CliqStatus
# Clique status message enumerated type with status.
#
# DevNotes
# - Temporary convert to Symbol to support #459 consolidation effort
# - Long term enum looks like a good idea (see #744)
# """
# @enum CliqStatus NULL INITIALIZED UPSOLVED MARGINALIZED DOWNSOLVED UPRECYCLED ERROR_STATUS
const CliqStatus = Symbol
## Currently same name by used as Symbol, e.g. :NULL, ...
## Older status names
# :null; :upsolved; :downsolved; :marginalized; :uprecycled,
## FIXME, consolidate at end of #459 work

# Used for UPWARD_DIFFERENTIAL, UPWARD_COMMON, DOWNWARD_COMMON marginalized types
abstract type MessagePassDirection end
struct UpwardPass <: MessagePassDirection end
struct DownwardPass <: MessagePassDirection end

abstract type MessageType end
struct NonparametricMessage <: MessageType end
struct ParametricMessage <: MessageType end

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
  softtype::T # rename to VariableType, DFG #603 
  # TODO -- DEPRECATE
  manifolds::Tuple{Vararg{Symbol}} # NOTE added during #459 effort
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

function compare(t1::TreeBelief, t2::TreeBelief)
  TP = true
  TP = TP && norm(t1.val - t2.val) < 1e-5
  TP = TP && norm(t1.bw - t2.bw) < 1e-5
  TP = TP && abs(t1.inferdim - t2.inferdim) < 1e-5
  TP = TP && t1.softtype == t2.softtype
  return TP
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
mutable struct LikelihoodMessage{T <: MessageType} <: AbstractPrior
  status::CliqStatus
  belief::Dict{Symbol, TreeBelief}
  variableOrder::Vector{Symbol}
  cliqueLikelihood::Union{Nothing,SamplableBelief}
  msgType::T
end


LikelihoodMessage(;status::CliqStatus=:NULL,
                   beliefDict::Dict=Dict{Symbol, TreeBelief}(),
                   variableOrder::Vector{Symbol}=Symbol[],
                   cliqueLikelihood=nothing,
                   msgType::T=NonparametricMessage() ) where {T <: MessageType} =
        LikelihoodMessage{T}(status, beliefDict, variableOrder, cliqueLikelihood, msgType)
#


function compare(l1::LikelihoodMessage,
                 l2::LikelihoodMessage;
                 skip::Vector{Symbol}=[] )
  #
  TP = true
  TP = TP && l1.status == l2.status
  TP = TP && l1.variableOrder == l2.variableOrder
  TP = TP && l1.msgType == l2.msgType
  TP = TP && l1.cliqueLikelihood |> typeof == l2.cliqueLikelihood |> typeof
  for (k,v) in l1.belief
    TP = TP && haskey(l2.belief, k)
    TP = TP && compare(v, l2.belief[k])
  end
end

==(l1::LikelihoodMessage,l2::LikelihoodMessage) = compare(l1,l2)


# FIXME, better standardize intermediate types
# used during nonparametric CK preparation, when information from multiple siblings must be shared together
const IntermediateSiblingMessages = Vector{Tuple{BallTreeDensity,Float64}}
const IntermediateMultiSiblingMessages = Dict{Symbol, IntermediateSiblingMessages}

const TempUpMsgPlotting = Dict{Symbol,Vector{Tuple{Symbol, Int, BallTreeDensity, Float64}}}


"""
$(TYPEDEF)
"""
mutable struct PotProd
    Xi::Symbol # Int
    prev::Array{Float64,2}
    product::Array{Float64,2}
    potentials::Array{BallTreeDensity,1}
    potentialfac::Vector{Symbol}
end
"""
$(TYPEDEF)
"""
mutable struct CliqGibbsMC
    prods::Array{PotProd,1}
    lbls::Vector{Symbol}
    CliqGibbsMC() = new()
    CliqGibbsMC(a,b) = new(a,b)
end
"""
$(TYPEDEF)

TO BE DEPRECATED
"""
mutable struct DebugCliqMCMC
  mcmc::Union{Nothing, Array{CliqGibbsMC,1}}
  outmsg::LikelihoodMessage
  outmsglbls::Dict{Symbol, Symbol} # Int
  priorprods::Vector{CliqGibbsMC}
  DebugCliqMCMC() = new()
  DebugCliqMCMC(a,b,c,d) = new(a,b,c,d)
end


"""
$(TYPEDEF)

TO BE DEPRECATED AND CONSOLIDATED
"""
mutable struct DownReturnBPType
  dwnMsg::LikelihoodMessage
  dbgDwn::DebugCliqMCMC
  IDvals::Dict{Symbol,TreeBelief}
  keepdwnmsgs::LikelihoodMessage
end






#
