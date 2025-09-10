
"""
$(TYPEDEF)

Message prior on all dimensions of a variable node in the factor graph.

Notes
- Only temporary existance during CSM operations.
"""
struct MsgPrior{T <: SamplableBelief} <: AbstractPriorObservation
  Z::T
  infoPerCoord::Vector{Float64}
  M::Any
end

# MsgPrior{T}() where {T} = new{T}()
# MsgPrior{T}(z::T, infd::R) where {T <: SamplableBelief, R <: Real} = new{T}(z, infd)
# function MsgPrior(z::T, infd::R) where {T <: SamplableBelief, R <: Real}
#     MsgPrior{T}(z, infd)
# end
function getSample(cf::CalcFactor{<:MsgPrior})
  return rand(cf.factor.Z, 1)
end

#TODO check these for manifolds, may need updating to samplePoint
# MKD already returns a vector of points
function getSample(cf::CalcFactor{<:MsgPrior{<:ManifoldKernelDensity}})
  mkd = cf.factor.Z
  return samplePoint(mkd.manifold, mkd)
end

getManifold(mp::MsgPrior{<:ManifoldKernelDensity}) = mp.Z.manifold
getManifold(mp::MsgPrior) = mp.M

#FIXME this will not work on manifolds
(cfo::CalcFactor{<:MsgPrior})(z, x1) = z .- x1

Base.@kwdef struct PackedMsgPrior <: AbstractPackedObservation
  Z::PackedBelief
  infoPerCoord::Vector{Float64}
end

function convert(::Type{PackedMsgPrior}, d::MsgPrior)
  return PackedMsgPrior(convert(PackedBelief, d.Z), d.infoPerCoord)
end
function convert(::Type{<:MsgPrior}, d::PackedMsgPrior)
  return MsgPrior(convert(SamplableBelief, d.Z), d.infoPerCoord)
end
