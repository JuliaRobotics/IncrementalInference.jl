
"""
$(TYPEDEF)

Message prior on all dimensions of a variable node in the factor graph.

Notes
- Only temporary existance during CSM operations.
"""
struct MsgPrior{T <: SamplableBelief} <: AbstractPrior
  Z::T
  infoPerCoord::Vector{Float64}
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


(cfo::CalcFactor{<:MsgPrior})(z, x1) = z .- x1



struct PackedMsgPrior <: AbstractPackedFactor
  Z::String
  infoPerCoord::Vector{Float64}
end

function convert(::Type{PackedMsgPrior}, d::MsgPrior)
  PackedMsgPrior(convert(String, d.Z), d.infoPerCoord) # TODO PackedSamplableBelief
end
function convert(::Type{<:MsgPrior}, d::PackedMsgPrior)
  MsgPrior(convert(SamplableBelief, d.Z), d.infoPerCoord)
end

