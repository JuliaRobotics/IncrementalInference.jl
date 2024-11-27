
## Distributions to JSON/Packed types

packDistribution(dtr::Categorical) = PackedCategorical(; p = dtr.p)
packDistribution(dtr::Uniform) = PackedUniform(; a = dtr.a, b = dtr.b)
packDistribution(dtr::Normal) = PackedNormal(; mu = dtr.μ, sigma = dtr.σ)
packDistribution(dtr::ZeroMeanDiagNormal) = PackedZeroMeanDiagNormal(; diag = dtr.Σ.diag)
packDistribution(dtr::ZeroMeanFullNormal) = PackedZeroMeanFullNormal(; cov = dtr.Σ.mat[:])
packDistribution(dtr::DiagNormal) = PackedDiagNormal(; mu = dtr.μ, diag = dtr.Σ.diag)
packDistribution(dtr::FullNormal) = PackedFullNormal(; mu = dtr.μ, cov = dtr.Σ.mat[:])
packDistribution(dtr::Rayleigh) = PackedRayleigh(; sigma = dtr.σ)

function packDistribution(dtr::AliasingScalarSampler)
  return PackedAliasingScalarSampler(; domain = dtr.domain, weights = dtr.weights.values)
end


## Unpack JSON/Packed to Distribution types

unpackDistribution(dtr::PackedCategorical) = Categorical(dtr.p ./ sum(dtr.p))
unpackDistribution(dtr::PackedUniform) = Uniform(dtr.a, dtr.b)
unpackDistribution(dtr::PackedNormal) = Normal(dtr.mu, dtr.sigma)
function unpackDistribution(dtr::PackedZeroMeanDiagNormal)
  return MvNormal(LinearAlgebra.Diagonal(map(abs2, sqrt.(dtr.diag))))
end # sqrt.(dtr.diag)
function unpackDistribution(dtr::PackedZeroMeanFullNormal)
  d = round(Int,sqrt(size(dtr.cov)[1]))
  return MvNormal(reshape(dtr.cov, d, d))
end
unpackDistribution(dtr::PackedDiagNormal) = MvNormal(dtr.mu, sqrt.(dtr.diag))
function unpackDistribution(dtr::PackedFullNormal)
  return MvNormal(dtr.mu, reshape(dtr.cov, length(dtr.mu), :))
end
unpackDistribution(dtr::PackedRayleigh) = Rayleigh(dtr.sigma)

function unpackDistribution(dtr::PackedAliasingScalarSampler)
  return AliasingScalarSampler(dtr.domain, dtr.weights ./ sum(dtr.weights))
end


# ## strip field from NamedTuple

# function _delete( nt::Union{<:NamedTuple, <:Dict{K,T}}, 
#                   key::K=:_type ) where {K,T}
#   #
#   kys = keys(nt)
#   # rm index
#   ridx = findfirst(k->k==key, kys)
#   # keep indices
#   idxs = setdiff(1:length(nt), ridx)
#   # to Dict
#   dict = OrderedDict{K,Any}()
#   for id in idxs
#     ky = kys[id]
#     dict[ky] = nt[ky]
#   end
#   # 

#   NamedTuple{Tuple(keys(dict))}(values(dict))
# end

## ===========================================================================================
## FIXME, should be obsolete and must be removed
## ===========================================================================================

# NOTE part of new effort to overhaul the SamplableBelief serialization approach
function convert(::Type{<:PackedSamplableBelief}, obj::StringThemSamplableBeliefs)
  return packDistribution(obj)
end
convert(::Type{<:SamplableBelief}, obj::PackedSamplableBelief) = unpackDistribution(obj)

function convert(::Type{<:PackedSamplableBelief}, nt::Union{NamedTuple, JSON3.Object})
  distrType = DFG.getTypeFromSerializationModule(nt._type)
  return distrType(; nt...)
end

##===================================================================================

# FIXME ON FIRE, must deprecate nested JSON written fields in all serialization
# TODO is string necessary, because unpacking templated e.g. PackedType{T} has problems, see DFG #668
function convert(::Type{String}, dtr::StringThemSamplableBeliefs)
  return JSON3.write(packDistribution(dtr))
end

function convert(::Type{<:SamplableBelief}, str_obj::AbstractString)
  #

  # go from stringified to generic packed (no type info)
  _pck = JSON3.read(str_obj)
  # NOTE, get the packed type from strong assumption that field `_type` exists in the 
  T = DFG.getTypeFromSerializationModule(_pck._type)
  # unpack again to described packedType
  pckT = JSON3.read(str_obj, T)

  # unpack to regular <:SamplableBelief
  return unpackDistribution(pckT)
end

#
