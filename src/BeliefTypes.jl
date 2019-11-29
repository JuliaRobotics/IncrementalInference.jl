
### SOME CONVERGENCE REQUIRED ---


"""
$(TYPEDEF)

Condensed representation of KernelDensityEstimate, by saving points and bandwidth
"""
mutable struct EasyMessage{T <: Tuple}
  pts::Array{Float64,2}
  bws::Array{Float64,1}
  manifolds::T
  inferdim::Float64
  EasyMessage{T}() where {T <: Tuple} = new{T}()
  EasyMessage{T}(a::Array{Float64,2}, b::Array{Float64,1}, manis::T, inferdim::Union{Float64, Int32, Int64}=0.0) where {T <: Tuple} = new{T}(a,b, manis, Float64(inferdim))
  EasyMessage{T}(p::BallTreeDensity, manis::T, inferdim::Union{Float64, Int32, Int64}=0.0) where {T <: Tuple}  = new{T}(getPoints(p), getBW(p)[:,1], manis, Float64(inferdim))
end
EasyMessage(a::Array{Float64,2}, b::Array{Float64,1}, manis::T, inferdim::Union{Float64, Int32, Int64}=0) where {T <: Tuple} = EasyMessage{T}(a, b, manis, inferdim)
EasyMessage(p::BallTreeDensity, manis::T, inferdim::Union{Float64, Int32, Int64}=0) where {T <: Tuple} = EasyMessage{T}(p, manis, inferdim)


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
  p::Dict{Symbol, EasyMessage}
end

### SOME CONVERGENCE REQUIRED ^^^

#DEV NOTE it looks like it can be consolidated into one type
# if we can pass messages similar to EasyMessage:
# pts::Array{Float64,2}
# bws::Array{Float64,1}
# option a tuple
# bellief::Dict{Symbol, NamedTuple{(:vec, :bw, :inferdim),Tuple{Array{Int64,1},Array{Int64,1},Float64}}}
# or an extra type
# or the MsgPrior/PackedMessagePrior, depending on the serialization requirement of the channel
# but I would think only one message type

# mutable struct NBPMessage <: Singleton
#   status::Symbol # Ek kort die in die boodskap
#   p::Dict{Symbol, EasyMessage}
# end

struct TreeBelief
  val::Array{Float64,2}
  bw::Array{Float64,2}
  inferdim::Float64
  # manifolds::T TODO ? JT not needed as you have the variable with all the info in it?
end
TreeBelief(p::BallTreeDensity, inferdim::Real=0.0) = TreeBelief(getPoints(p), getBW(p), inferdim)
TreeBelief(val::Array{Float64,2}, bw::Array{Float64,2}, inferdim::Real=0.0) = TreeBelief(val, bw, inferdim)

abstract type AbstractBeliefMessage end
#JT Ek maak nog whahahaha, maar kom ons inheret van AbstractBeliefMessage?
# of ons gebruik T parameter in belief::Dict{Symbol, T} where T <: Union{BallTreeDensity, Vector{Float64}}
struct ParametricBeliefMessage <: AbstractBeliefMessage
  # TODO JT nog velde, maar dis wat ek sover aan kan dink
  status::Symbol
  belief::Dict{Symbol, TreeBelief}
  # of dalk named tuple maak
  # bellief::Dict{Symbol, NamedTuple{(:val, :bw, :inferdim),Tuple{Array{Int64,1},Array{Int64,1},Float64}}}
end

ParametricBeliefMessage(status::Symbol) =
        ParametricBeliefMessage(status, Dict{Symbol, TreeBelief}())
