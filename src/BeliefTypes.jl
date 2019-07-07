
### SOME CONVERGENCE REQUIRED ---


"""
$(TYPEDEF)

Condensed representation of KernelDensityEstimate, by saving points and bandwidth
"""
mutable struct EasyMessage{T <: Tuple}
  pts::Array{Float64,2}
  bws::Array{Float64,1}
  manifolds::T
  fulldim::Bool
  EasyMessage{T}() where {T <: Tuple} = new{T}()
  EasyMessage{T}(a::Array{Float64,2}, b::Array{Float64,1}, manis::T, fulldim::Bool=true) where {T <: Tuple} = new{T}(a,b, manis, fulldim)
  EasyMessage{T}(p::BallTreeDensity, manis::T, fulldim::Bool=true) where {T <: Tuple}  = new{T}(getPoints(p), getBW(p)[:,1], manis, fulldim)
end
EasyMessage(a::Array{Float64,2}, b::Array{Float64,1}, manis::T, fulldim::Bool=true) where {T <: Tuple} = EasyMessage{T}(a, b, manis, fulldim)
EasyMessage(p::BallTreeDensity, manis::T, fulldim::Bool=true) where {T <: Tuple} = EasyMessage{T}(p, manis, fulldim)

const TempBeliefMsg = Dict{Symbol, Tuple{BallTreeDensity, Vector{Bool}}}
"""
$(TYPEDEF)
"""
mutable struct NBPMessage <: Singleton
  p::Dict{Symbol, EasyMessage}
end

### SOME CONVERGENCE REQUIRED ^^^
