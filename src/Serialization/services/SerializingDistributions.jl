
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

# A more specialized constructor to help serialization processes reading Any[Any[1,2,3,..]] rather than floats. 
function PackedHeatmapGridDensity(
  _type::String,
  data::AbstractVector, # {Any} 
  domain::AbstractVector, # {Any}
  hint_callback::String,
  bw_factor::Float64,
  N::Int64,
)
  #
  # TODO data might not be type Float64, should store and recover as performance enhancement (if user specified different element type)
  data_ = Vector{Vector{Float64}}(undef, length(data))
  for (i, dat) in enumerate(data)
    dat_ = replace(dat, nothing => 0)
    data_[i] = float.(dat_)
  end
  domain_ = tuple(float.(domain[1]), float.(domain[2]))

  return PackedHeatmapGridDensity(_type, data_, domain_, hint_callback, bw_factor, N)
end

function packDistribution(obj::HeatmapGridDensity)
  #
  data_ = obj.data
  @cast data[j][i] := data_[i, j]
  # str = convert(SamplableBelief, obj.densityFnc)
  N = Npts(obj.densityFnc)
  # TODO misses the hint...
  return PackedHeatmapGridDensity(
    "IncrementalInference.PackedHeatmapGridDensity",
    data,
    obj.domain,
    "",
    obj.bw_factor,
    N,
  )
end

function packDistribution(dtr::LevelSetGridNormal)
  return PackedLevelSetGridNormal(
    "IncrementalInference.PackedLevelSetGridNormal",
    dtr.level,
    dtr.sigma,
    dtr.sigma_scale,
    convert(PackedHeatmapGridDensity, dtr.heatmap),
  )
end
#

function parchDistribution(hgd::HeatmapGridDensity)
  @assert 2 <= size(hgd.data, 1) "parchDistribution of HeatmapGridDensity can only be done when `.data` is larger than 2x1"

  data = Matrix{eltype(hgd.data)}(undef, 2, 2)
  data[1, 1] = hgd.data[1, 1]
  # data[2,2] = hgd.data[2,2] # disable since data might be a single column in unusual cases
  data[2, 1] = size(hgd.data, 1)
  data[1, 2] = size(hgd.data, 2)

  domain = hgd.domain
  hint_callback = hgd.hint_callback
  bw_factor = hgd.bw_factor
  densityFnc = parchDistribution(hgd.densityFnc)

  return HeatmapGridDensity(data, domain, hint_callback, bw_factor, densityFnc)
end

## Unpack JSON/Packed to Distribution types

unpackDistribution(dtr::PackedCategorical) = Categorical(dtr.p ./ sum(dtr.p))
unpackDistribution(dtr::PackedUniform) = Uniform(dtr.a, dtr.b)
unpackDistribution(dtr::PackedNormal) = Normal(dtr.mu, dtr.sigma)
function unpackDistribution(dtr::PackedZeroMeanDiagNormal)
  return MvNormal(LinearAlgebra.Diagonal(map(abs2, sqrt.(dtr.diag))))
end # sqrt.(dtr.diag)
function unpackDistribution(dtr::PackedZeroMeanFullNormal)
  return MvNormal(reshape(dtr.cov, length(dtr.mu), :))
end
unpackDistribution(dtr::PackedDiagNormal) = MvNormal(dtr.mu, sqrt.(dtr.diag))
function unpackDistribution(dtr::PackedFullNormal)
  return MvNormal(dtr.mu, reshape(dtr.cov, length(dtr.mu), :))
end
unpackDistribution(dtr::PackedRayleigh) = Rayleigh(dtr.sigma)

function unpackDistribution(dtr::PackedAliasingScalarSampler)
  return AliasingScalarSampler(dtr.domain, dtr.weights ./ sum(dtr.weights))
end

function unpackDistribution(obj::PackedHeatmapGridDensity)
  #
  # do intermediate conversions
  data_ = obj.data
  data__ = map(x -> collect(x), data_)
  @cast data[i, j] := data__[j][i]
  _data__ = collect(data)
  # densFnc = convert(SamplableBelief, obj.densityFnc)
  # build the final object, misses the hint...
  return HeatmapGridDensity(
    _data__,
    obj.domain,
    obj.hint_callback == "" ? nothing : nothing,
    obj.bw_factor;
    N = obj.N,
  )
end

function unpackDistribution(dtr::PackedLevelSetGridNormal)
  return LevelSetGridNormal(
    dtr.level,
    dtr.sigma,
    dtr.sigma_scale,
    convert(HeatmapGridDensity, dtr.heatmap),
  )
end
#

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
## Converts still necessary?
## ===========================================================================================

# NOTE part of new effort to overhaul the SamplableBelief serialization approach
function convert(::Type{<:PackedSamplableBelief}, obj::StringThemSamplableBeliefs)
  return packDistribution(obj)
end
convert(::Type{<:SamplableBelief}, obj::PackedSamplableBelief) = unpackDistribution(obj)

function convert(::Type{<:PackedSamplableBelief}, nt::NamedTuple)
  distrType = DFG.getTypeFromSerializationModule(nt._type)
  return distrType(; nt...)
end

##===================================================================================

# FIXME ON FIRE, must deprecate nested JSON written fields in all serialization
# TODO is string necessary, because unpacking templated e.g. PackedType{T} has problems, see DFG #668
function convert(::Type{String}, dtr::StringThemSamplableBeliefs)
  return JSON2.write(packDistribution(dtr))
end

function convert(::Type{<:SamplableBelief}, str_obj::AbstractString)
  #

  # go from stringified to generic packed (no type info)
  _pck = JSON2.read(str_obj)
  # NOTE, get the packed type from strong assumption that field `_type` exists in the 
  T = DFG.getTypeFromSerializationModule(_pck[:_type])
  # unpack again to described packedType
  pckT = JSON2.read(str_obj, T)

  # unpack to regular <:SamplableBelief
  return unpackDistribution(pckT)
end

#
