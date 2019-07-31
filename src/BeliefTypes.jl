
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
"""
mutable struct NBPMessage <: Singleton
  p::Dict{Symbol, EasyMessage}
end

### SOME CONVERGENCE REQUIRED ^^^
