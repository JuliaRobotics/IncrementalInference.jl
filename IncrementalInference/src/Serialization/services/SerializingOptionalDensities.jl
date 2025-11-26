



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

function DFG.pack(obj::HeatmapGridDensity)
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

function DFG.pack(dtr::LevelSetGridNormal)
  return PackedLevelSetGridNormal(
    "IncrementalInference.PackedLevelSetGridNormal",
    dtr.level,
    dtr.sigma,
    dtr.sigma_scale,
    DFG.pack(dtr.heatmap),
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

function DFG.unpack(obj::PackedHeatmapGridDensity)
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

function DFG.unpack(dtr::PackedLevelSetGridNormal)
  return LevelSetGridNormal(
    dtr.level,
    dtr.sigma,
    dtr.sigma_scale,
    unpack(dtr.heatmap),
  )
end


#