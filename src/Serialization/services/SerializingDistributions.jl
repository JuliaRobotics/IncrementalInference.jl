

## Distributions to JSON/Packed types

packDistribution(dtr::Categorical) = PackedCategorical(; p=dtr.p )
packDistribution(dtr::Uniform) = PackedUniform(; a=dtr.a, b=dtr.b )
packDistribution(dtr::Normal) = PackedNormal(; mu=dtr.μ, sigma=dtr.σ )
packDistribution(dtr::ZeroMeanDiagNormal) = PackedZeroMeanDiagNormal(; diag=dtr.Σ.diag )
packDistribution(dtr::ZeroMeanFullNormal) = PackedZeroMeanFullNormal(; cov=dtr.Σ.mat[:] )
packDistribution(dtr::DiagNormal) = PackedDiagNormal(; mu=dtr.μ, diag=dtr.Σ.diag )
packDistribution(dtr::FullNormal) = PackedFullNormal(; mu=dtr.μ, cov=dtr.Σ.mat[:] )
packDistribution(dtr::Rayleigh) = PackedRayleigh(; sigma=dtr.σ )

packDistribution(dtr::AliasingScalarSampler) = PackedAliasingScalarSampler(; domain=dtr.domain, weights=dtr.weights.values )

# A more specialized constructor to help serialization processes reading Any[Any[1,2,3,..]] rather than floats. 
function PackedHeatmapGridDensity(
    _type::String, 
    data::AbstractVector, # {Any} 
    domain::AbstractVector, # {Any}
    hint_callback::String, 
    bw_factor::Float64, 
    N::Int64 )
  #
  data_ = Vector{Vector{Float64}}(undef, length(data))
  for (i,dat) in enumerate(data)
    data_[i] = float.(dat)
  end
  domain_ = tuple(float.(domain[1]), float.(domain[2]))

  PackedHeatmapGridDensity(_type, data_, domain_, hint_callback, bw_factor, N)
end

function packDistribution( obj::HeatmapGridDensity )
  #
  data_ = obj.data
  @cast data[j][i] := data_[i,j]
  # str = convert(SamplableBelief, obj.densityFnc)
  N = Npts(obj.densityFnc)
  # TODO misses the hint...
  PackedHeatmapGridDensity( "IncrementalInference.PackedHeatmapGridDensity",
                            data,
                            obj.domain,
                            "", 
                            obj.bw_factor,
                            N )
end

packDistribution(dtr::LevelSetGridNormal) = PackedLevelSetGridNormal( "IncrementalInference.PackedLevelSetGridNormal",
                                                                      dtr.level,
                                                                      dtr.sigma,
                                                                      dtr.sigma_scale,
                                                                      convert(PackedHeatmapGridDensity, dtr.heatmap) )
#


## Unpack JSON/Packed to Distribution types

unpackDistribution(dtr::PackedCategorical) = Categorical( dtr.p ./ sum(dtr.p) )
unpackDistribution(dtr::PackedUniform) = Uniform(dtr.a, dtr.b )
unpackDistribution(dtr::PackedNormal) = Normal( dtr.mu, dtr.sigma )
unpackDistribution(dtr::PackedZeroMeanDiagNormal) = MvNormal( sqrt.(dtr.diag) ) # LinearAlgebra.Diagonal(map(abs2, σ))
unpackDistribution(dtr::PackedZeroMeanFullNormal) = MvNormal( reshape(dtr.cov, length(dtr.mu), :) )
unpackDistribution(dtr::PackedDiagNormal) = MvNormal( dtr.mu, sqrt.(dtr.diag) )
unpackDistribution(dtr::PackedFullNormal) = MvNormal( dtr.mu, reshape(dtr.cov, length(dtr.mu), :) )
unpackDistribution(dtr::PackedRayleigh) = Rayleigh( dtr.sigma )

unpackDistribution(dtr::PackedAliasingScalarSampler) = AliasingScalarSampler( dtr.domain, dtr.weights ./ sum(dtr.weights) )

function unpackDistribution(obj::PackedHeatmapGridDensity)
  #
  # do intermediate conversions
  data_ = obj.data
  data__ = map(x->collect(x), data_)
  @cast data[i,j] := data__[j][i]
  _data__ = collect(data)
  # densFnc = convert(SamplableBelief, obj.densityFnc)
  # build the final object, misses the hint...
  HeatmapGridDensity( _data__,
                      obj.domain,
                      obj.hint_callback == "" ? nothing : nothing,
                      obj.bw_factor;
                      N=obj.N )
end

unpackDistribution(dtr::PackedLevelSetGridNormal) = LevelSetGridNormal( dtr.level,
                                                                        dtr.sigma,
                                                                        dtr.sigma_scale,
                                                                        convert(HeatmapGridDensity, dtr.heatmap) )
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
convert(::Type{<:PackedSamplableBelief}, obj::StringThemSamplableBeliefs) = packDistribution(obj)
convert(::Type{<:SamplableBelief}, obj::PackedSamplableBelief) = unpackDistribution(obj)

function convert(::Type{<:PackedSamplableBelief}, nt::NamedTuple)
  distrType = DFG.getTypeFromSerializationModule(nt._type)
  distrType(;nt...)
end

##===================================================================================


# FIXME ON FIRE, must deprecate nested JSON written fields in all serialization
# TODO is string necessary, because unpacking templated e.g. PackedType{T} has problems, see DFG #668
convert(::Type{String}, dtr::StringThemSamplableBeliefs) = JSON2.write(packDistribution(dtr))

function convert(::Type{<:SamplableBelief}, str_obj::AbstractString)
  #

  # go from stringified to generic packed (no type info)
  _pck = JSON2.read(str_obj)
  # NOTE, get the packed type from strong assumption that field `_type` exists in the 
  T = DFG.getTypeFromSerializationModule( _pck[:_type] )  
  # unpack again to described packedType
  pckT = JSON2.read(str_obj, T)

  # unpack to regular <:SamplableBelief
  return unpackDistribution(pckT)
end


#