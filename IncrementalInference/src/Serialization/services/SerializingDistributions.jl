
function pack(dtr::AliasingScalarSampler)
  return PackedAliasingScalarSampler(; domain = dtr.domain, weights = dtr.weights.values)
end

function unpack(dtr::PackedAliasingScalarSampler)
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
function convert(::Type{<:PackedBelief}, obj::StringThemSamplableBeliefs)
  return pack(obj)
end
convert(::Type{<:SamplableBelief}, obj::PackedBelief) = unpack(obj)

##===================================================================================

# FIXME ON FIRE, must deprecate nested JSON written fields in all serialization
# TODO is string necessary, because unpacking templated e.g. PackedType{T} has problems, see DFG #668
function convert(::Type{String}, dtr::StringThemSamplableBeliefs)
  return JSON3.write(pack(dtr))
end

function convert(::Type{<:SamplableBelief}, str_obj::AbstractString)
  #

  # go from stringified to generic packed (no type info)
  _pck = JSON3.read(str_obj)
  # NOTE, get the packed type from strong assumption that field `_type` exists in the 
  T = getTypeFromSerializationModule(_pck._type)
  # unpack again to described packedType
  pckT = JSON3.read(str_obj, T)

  # unpack to regular <:SamplableBelief
  return unpack(pckT)
end

#
