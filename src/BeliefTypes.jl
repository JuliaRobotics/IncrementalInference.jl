
export TreeBelief, LikelihoodMessage

# const TreeBelief{T} = EasyMessage{T}
# const EasyMessage{T} = TreeBelief{T}

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
    CliqStatus
Clique status message enumerated type with status:
initialized, upsolved, marginalized, downsolved, uprecycled
"""
@enum CliqStatus initialized upsolved marginalized downsolved uprecycled error_status


"""
  $(TYPEDEF)
Belief message for message passing on the tree.

Notes:
- belief -> common mode
- cobelief -> differential mode

DevNotes:
- Objective: `MvNormal(μ=[:x0;:x2;:l5], Σ=[+ * *; * + *; * * +])`

  $(TYPEDFIELDS)
"""
struct LikelihoodMessage
  status::CliqStatus
  belief::Dict{Symbol, TreeBelief}
  cobelief::NamedTuple{(:varlbl, :μ, :Σ),Tuple{Vector{Symbol}, Vector{Float64}, Matrix{Float64}}} #TODO name something mathier
end

LikelihoodMessage(status::CliqStatus) =
        LikelihoodMessage(status, Dict{Symbol, TreeBelief}(), (varlbl=Symbol[], μ=Float64[], Σ=Matrix{Float64}(undef,0,0)))

LikelihoodMessage(status::CliqStatus, cobelief) =
        LikelihoodMessage(status, Dict{Symbol, TreeBelief}(), cobelief)


#

const BeliefMessage = LikelihoodMessage




### EVERYTHING BELOW IS/SHOULD BE DEPRECATED



"""
$(TYPEDEF)

Condensed representation of KernelDensityEstimate, by saving points and bandwidth
"""
mutable struct EasyMessage{T <: Tuple}
  val::Array{Float64,2}
  bw::Array{Float64,1}
  manifolds::T
  inferdim::Float64
  EasyMessage{T}() where {T <: Tuple} = new{T}()
  EasyMessage{T}(a::Array{Float64,2}, b::Array{Float64,1}, manis::T, inferdim::Union{Float64, Int32, Int64}=0.0) where {T <: Tuple} = new{T}(a,b, manis, Float64(inferdim))
  EasyMessage{T}(p::BallTreeDensity, manis::T, inferdim::Union{Float64, Int32, Int64}=0.0) where {T <: Tuple}  = new{T}(getPoints(p), getBW(p)[:,1], manis, Float64(inferdim))
end

EasyMessage(a::Array{Float64,2}, b::Array{Float64,1}, manis::T, inferdim::Union{Float64, Int32, Int64}=0) where {T <: Tuple} = EasyMessage{T}(a, b, manis, inferdim)
EasyMessage(p::BallTreeDensity, manis::T, inferdim::Union{Float64, Int32, Int64}=0) where {T <: Tuple} = EasyMessage{T}(p, manis, inferdim)




# Deprecated, replaced by LikelihoodMessage
# TODO - remove
const TempBeliefMsg = Dict{Symbol, Tuple{BallTreeDensity, Float64}}

# Dict{Symbol,   -- is for variable label
#  Vector{       -- multiple msgs for the same variable
#   Symbol,      -- Clique index
#   Int,         -- Depth in tree
#   BTD          -- Belief estimate
#   inferredDim  -- Information count
#  }
const TempUpMsgPlotting = Dict{Symbol,Vector{Tuple{Symbol, Int, BallTreeDensity, Float64}}}



"""
$(TYPEDEF)

DESPARATELY NEEDS TO BE UPDATED TO USE TempBeliefMsg DEFINITION (start of refactor).
"""
mutable struct NBPMessage <: Singleton
  p::Union{Dict{Symbol, TreeBelief}, Dict{Symbol, EasyMessage}}
end




#
