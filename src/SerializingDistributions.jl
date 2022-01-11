
## Pack Distributions to JSON/Packed types

packDistribution(dtr::Categorical) = PackedCategorical(; p=dtr.p )
packDistribution(dtr::Uniform) = PackedUniform(; a=dtr.a, b=dtr.b )
packDistribution(dtr::Normal) = PackedNormal(; mu=dtr.μ, sigma=dtr.σ )
packDistribution(dtr::ZeroMeanDiagNormal) = PackedZeroMeanDiagNormal(; diag=dtr.Σ.diag )
packDistribution(dtr::ZeroMeanFullNormal) = PackedZeroMeanFullNormal(; cov=dtr.Σ.mat )
packDistribution(dtr::DiagNormal) = PackedDiagNormal(; mu=dtr.μ, diag=dtr.Σ.diag )
packDistribution(dtr::FullNormal) = PackedFullNormal(; mu=dtr.μ, cov=dtr.Σ.mat )

packDistribution(dtr::AliasingScalarSampler) = PackedAliasingScalarSampler(; domain=dtr.domain, weights=dtr.weights.values )

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


## Unpack distribution JSON/Packed types

unpackDistribution(dtr::PackedCategorical) = Categorical( dtr.p ./ sum(dtr.p) )
unpackDistribution(dtr::PackedUniform) = Uniform(dtr.a, dtr.b )
unpackDistribution(dtr::PackedNormal) = Normal( dtr.mu, dtr.sigma )
unpackDistribution(dtr::PackedZeroMeanDiagNormal) = MvNormal( sqrt.(dtr.diag) )
unpackDistribution(dtr::PackedZeroMeanFullNormal) = MvNormal( dtr.cov )
unpackDistribution(dtr::PackedDiagNormal) = MvNormal( dtr.mu, sqrt.(dtr.diag) )
unpackDistribution(dtr::PackedFullNormal) = MvNormal( dtr.mu, dtr.cov )

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



## ===========================================================================================
## Converts still necessary?
## ===========================================================================================


# function convert(::Type{<:PackedSamplableBelief}, obj::Distributions.Distribution)
  
# end


function convert( ::Union{Type{<:PackedSamplableBelief},Type{<:PackedUniform}},
                  obj::Distributions.Uniform )
  #
  packed = packDistribution(obj)

  # FIXME remove JSON writing here! 
  return JSON2.write(packed)
end

convert(::Type{<:SamplableBelief}, obj::PackedUniform) = unpackDistribution(obj)


# NOTE part of new effort to overhaul the SamplableBelief serialization approach
# maybe it all becomes a JSON struct sort of thing in the long run.
# convert(::Type{<:PackedSamplableBelief}, obj::StringThemSamplableBeliefs) = string(obj)


# FIXME DEPRECATE TO BETTER JSON with ._type field STANDARD
function convert(::Type{<:PackedSamplableBelief}, obj::SamplableBelief)
  # FIXME, prep for switch
  packDistribution(obj)
  
  # FIXME must use string, because unpacking templated e.g. PackedType{T} has problems, see DFG #668
  string(obj)
end


convert(::Union{Type{<:SamplableBelief},Type{<:HeatmapGridDensity}},
        obj::PackedHeatmapGridDensity) = unpackDistribution(obj)



convert(::Union{Type{<:PackedSamplableBelief},Type{<:PackedHeatmapGridDensity}}, 
        obj::HeatmapGridDensity ) = packDistribution(obj)
#

convert(::Union{Type{<:SamplableBelief},Type{<:LevelSetGridNormal}}, 
        obj::PackedLevelSetGridNormal) = unpackDistribution(obj)


convert(::Union{Type{<:PackedSamplableBelief},Type{<:PackedLevelSetGridNormal}}, 
        obj::LevelSetGridNormal) = packDistribution(obj)




# New features towards standardizing distribution serialization
# # Assumes DFG/IIF serialized distributions have a `PackedType._type::String = "MyModule.MyPackedDistributionDensityType"`
# # also see DFG #590
# function convert( ::Type{String}, 
#                   obj::PackedSamplableBelief )
#   #
#   _typ = DFG.getTypeFromSerializationModule(obj._type)
# end

#