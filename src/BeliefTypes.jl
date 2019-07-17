
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
  EasyMessage{T}(a::Array{Float64,2}, b::Array{Float64,1}, manis::T, inferdim::Float64=0) where {T <: Tuple} = new{T}(a,b, manis, inferdim)
  EasyMessage{T}(p::BallTreeDensity, manis::T, inferdim::Float64=0) where {T <: Tuple}  = new{T}(getPoints(p), getBW(p)[:,1], manis, inferdim)
end
EasyMessage(a::Array{Float64,2}, b::Array{Float64,1}, manis::T, inferdim::Float64=0) where {T <: Tuple} = EasyMessage{T}(a, b, manis, inferdim)
EasyMessage(p::BallTreeDensity, manis::T, inferdim::Float64=0) where {T <: Tuple} = EasyMessage{T}(p, manis, inferdim)

const TempBeliefMsg = Dict{Symbol, Tuple{BallTreeDensity, Float64}} # Dict{Symbol, Tuple{BallTreeDensity, Vector{Bool}}}
"""
$(TYPEDEF)
"""
mutable struct NBPMessage <: Singleton
  p::Dict{Symbol, EasyMessage}
end

### SOME CONVERGENCE REQUIRED ^^^
