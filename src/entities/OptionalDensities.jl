# entities immediately available as private members in IIF, but requires other packages for actual use

# only export once the convenience constructors are available along with conditional Interpolations dependency

struct FluxModelsDistribution{ID,OD,P,D<:AbstractArray}
  # shape of the input data
  inputDim::NTuple{ID,Int}
  # shape of the output data
  outputDim::NTuple{OD,Int}
  # actual Flux models
  models::Vector{P}
  # the data used for prediction, must be <: AbstractArray
  data::D
  # shuffle model predictions relative to particle index at each sampling
  shuffle::Base.RefValue{Bool}
  # false for default serialization with model info, set true for separate storage of models 
  serializeHollow::Base.RefValue{Bool}
    # # TODO remove requirement and standardize sampler API
    # specialSampler::Function
end



"""
    $TYPEDEF

Generate a `<:SamplableBelief` from a heatmap, e.g. a digital elevation model.

Notes
- Give in heatmap and grid, and object becomes a density function that can also be sampled.
- Sampling can be more nuanced by injecting a hint, or location of interest:
  - Mostly aimed at limiting compute when faced with massive heatmaps, e.g. nav units are 10's but map is ~1e6.
- Density approximation is constructed on Guassian measurement assumption of level set and sigma variation.
- Assume data is on a regular grid on TranslationGroup(2)

DevNotes:
- Generalize to scalar fields on any Manifold.
- Generalize to vector fields if interpolation is sensible.
"""
struct HeatmapDensityRegular{T <: Real, H <: Union{<:Function, Nothing}, B <: Union{ManifoldKernelDensity, BallTreeDensity}}
  """intensity data, on regular grid"""
  data::Matrix{T}
  """domain as grid or locations at which scalar intensity elements exist"""
  domain::Tuple{<:AbstractVector{T},<:AbstractVector{T}}
  """use location hint to focus sampling to specific area of data, requires additional info at `getSample`
      assumed the callback will return _____ NOT ACTIVE YET"""
  hint_callback::H
  """level at which to extract the set of interest"""
  level::T
  """one sigma value associated with measurement noise of `level` against `data`"""
  sigma::T
  """make samplible region of interest from data be `sigma_scale` from `level`, e.g. 3*sigma."""
  sigma_scale::T
  """general rule for kernel bandwidths used in construction of density, e.g. 0.7 of domain grid spacing"""
  bw_factor::T 
  """density function as samplable representation of the data over the domain"""
  densityFnc::B # TODO change to ::ManifoldKernelDensity{TranslationGroup(2),BallTreeDensity}
end

(hmd::HeatmapDensityRegular)(w...;kw...) = hmd.densityFnc(w...;kw...)

function sampleTangent(M::AbstractManifold, hms::HeatmapDensityRegular)
  sampleTangent(M, hms.densityFnc)
end


function Base.show(io::IO, x::HeatmapDensityRegular{T,H,B}) where {T,H,B}
  printstyled(io, "HeatmapDensityRegular{", bold=true, color=:blue)
  println(io)
  printstyled(io, "    T", color=:magenta, bold=true )
  println(io, "      = ", T)
  printstyled(io, "    H", color=:magenta, bold=true )
  println(io, "`int  = ", H )
  printstyled(io, "    B", color=:magenta, bold=true )
  println(io, "      = ", B )
  printstyled(io, " }", color=:blue, bold=true)
  println(io, "(")
  println(io, "  data:       ", size(x.data))
  println(io, "    min/max:    ", round(minimum(x.data),digits=5), " / ", round(maximum(x.data),digits=5))
  println(io, "  domain:     ", size(x.domain[1]), ", ", size(x.domain[2]))
  println(io, "    min/max:    ", round(minimum(x.domain[1]),digits=5), " / ", round(maximum(x.domain[1]),digits=5))
  println(io, "    min/max:    ", round(minimum(x.domain[2]),digits=5), " / ", round(maximum(x.domain[2]),digits=5))
  println(io, "  level:      ", x.level)
  println(io, "  sigma:      ", x.sigma)
  println(io, "  sig.scale:  ", x.sigma_scale)
  println(io, "  bw_factor:  ", x.bw_factor)
  print(io, "  ")
  show(io, x.densityFnc)
  nothing
end

Base.show(io::IO, ::MIME"text/plain", x::HeatmapDensityRegular) = show(io, x)

Base.show(io::IO, ::MIME"application/prs.juno.inline", x::HeatmapDensityRegular) = show(io,x)






#