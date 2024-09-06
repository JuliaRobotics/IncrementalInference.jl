# heatmap sampler (experimental)


(hmd::HeatmapGridDensity)(w...; kw...) = hmd.densityFnc(w...; kw...)

function sampleTangent(M::AbstractManifold, hms::HeatmapGridDensity)
  return sampleTangent(M, hms.densityFnc)
end

function Base.show(io::IO, x::HeatmapGridDensity{T, H, B}) where {T, H, B}
  printstyled(io, "HeatmapGridDensity{"; bold = true, color = :blue)
  println(io)
  printstyled(io, "    T"; color = :magenta, bold = true)
  println(io, "      = ", T)
  printstyled(io, "    H"; color = :magenta, bold = true)
  println(io, "`int  = ", H)
  printstyled(io, "    B"; color = :magenta, bold = true)
  println(io, "      = ", B)
  printstyled(io, " }"; color = :blue, bold = true)
  println(io, "(")
  println(io, "  data:       ", size(x.data))
  println(
    io,
    "    min/max:    ",
    round(minimum(x.data); digits = 5),
    " / ",
    round(maximum(x.data); digits = 5),
  )
  println(io, "  domain:     ", size(x.domain[1]), ", ", size(x.domain[2]))
  println(
    io,
    "    min/max:    ",
    round(minimum(x.domain[1]); digits = 5),
    " / ",
    round(maximum(x.domain[1]); digits = 5),
  )
  println(
    io,
    "    min/max:    ",
    round(minimum(x.domain[2]); digits = 5),
    " / ",
    round(maximum(x.domain[2]); digits = 5),
  )
  println(io, "  bw_factor:  ", x.bw_factor)
  print(io, "  ")
  show(io, x.densityFnc)
  return nothing
end

Base.show(io::IO, ::MIME"text/plain", x::HeatmapGridDensity) = show(io, x)
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::HeatmapGridDensity) = show(io, x)

"""
    $SIGNATURES

Internal function for updating HGD.  
  
Notes
- Likely to be used for [unstashing packed factors](@ref section_stash_unstash) via [`preambleCache`](@ref).
- Counterpart to `AMP._update!` function for stashing of either MKD or HGD.
"""
function _update!(
  dst::HeatmapGridDensity{T, H, B},
  src::HeatmapGridDensity{T, H, B},
) where {T, H, B}
  @assert size(dst.data) == size(src.data) "Updating HeatmapDensityGrid can only be done for data of the same size"
  dst.data .= src.data
  if !isapprox(dst.domain[1], src.domain[1])
    dst.domain[1] .= src.domain[1]
  end
  if !isapprox(dst.domain[2], src.domain[2])
    dst.domain[2] .= src.domain[2]
  end
  AMP._update!(dst.densityFnc, src.densityFnc)
  return dst
end


##

(lsg::LevelSetGridNormal)(w...; kw...) = lsg.densityFnc(w...; kw...)

function sampleTangent(M::AbstractManifold, lsg::LevelSetGridNormal)
  return sampleTangent(M, lsg.heatmap.densityFnc)
end

function Base.show(io::IO, x::LevelSetGridNormal{T, H}) where {T, H}
  printstyled(io, "LevelSetGridNormal{"; bold = true, color = :blue)
  println(io)
  printstyled(io, "    T"; color = :magenta, bold = true)
  println(io, "      = ", T)
  printstyled(io, "    H"; color = :magenta, bold = true)
  println(io, "`int  = ", H)
  printstyled(io, " }"; color = :blue, bold = true)
  println(io, "(")
  println(io, "  level:      ", x.level)
  println(io, "  sigma:      ", x.sigma)
  println(io, "  sig.scale:  ", x.sigma_scale)
  println(io, "  heatmap:    ")
  show(io, x.heatmap)
  return nothing
end

Base.show(io::IO, ::MIME"text/plain", x::LevelSetGridNormal) = show(io, x)
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::LevelSetGridNormal) = show(io, x)


##

getManifold(hgd::HeatmapGridDensity) = getManifold(hgd.densityFnc)
getManifold(lsg::LevelSetGridNormal) = getManifold(lsg.heatmap)

AMP.sample(hgd::HeatmapGridDensity, w...; kw...) = sample(hgd.densityFnc, w...; kw...)

"""
    $SIGNATURES

Get the grid positions at the specified height (within the provided spreads)

DevNotes
- TODO Should this be consolidated with AliasingScalarSampler? See IIF #1341
"""
function sampleHeatmap(
  roi::AbstractMatrix{<:Real},
  x_grid::AbstractVector{<:Real},
  y_grid::AbstractVector{<:Real},
  thres::Real = 1e-14,
)
  #

  # mask the region of interest above the sampling threshold value
  mask = thres .< roi

  idx2d = findall(mask)  # 2D indices
  pos = (v -> [x_grid[v[1]], y_grid[v[2]]]).(idx2d)
  weights = (v -> roi[v[1], v[2]]).(idx2d)
  weights ./= sum(weights)

  return pos, weights
end

# TODO make n-dimensional, and later on-manifold
# TODO better standardize for heatmaps on manifolds w MKD
function fitKDE(
  support,
  weights,
  x_grid::AbstractVector{<:Real},
  y_grid::AbstractVector{<:Real};
  bw_factor::Real = 0.7,
)
  #
  # 1. set the bandwidth 
  x_spacing = Statistics.mean(diff(x_grid))
  y_spacing = Statistics.mean(diff(y_grid))
  kernel_ = bw_factor * 0.5 * (x_spacing + y_spacing) # 70% of the average spacing
  kernel_bw = [kernel_; kernel_]                  # same bw in x and y
  # fit KDE
  return kde!(support, kernel_bw, weights)
end

# Helper function to construct HGD
function HeatmapGridDensity(
  # NAME SUGGESTION: SAMPLEABLE_FIELD, GenericField
  field_on_grid::AbstractMatrix{<:Real},
  domain::Tuple{<:AbstractVector{<:Real}, <:AbstractVector{<:Real}},
  # encapsulate above
  hint_callback::Union{<:Function, Nothing} = nothing,
  bw_factor::Real = 0.7;  # kde spread between domain points 
  N::Int = 10000,
)
  #
  pos, weights_ = sampleHeatmap(field_on_grid, domain..., 0)
  # recast to the appropriate shape
  @cast support_[i, j] := pos[j][i]

  # constuct a pre-density from which to draw intermediate samples
  # TODO remove extraneous collect()
  density_ = fitKDE(collect(support_), weights_, domain...; bw_factor = bw_factor)
  pts_preIS, = sample(density_, N)

  @cast vec_preIS[j][i] := pts_preIS[i, j]

  # weight the intermediate samples according to interpolation of raw field_on_grid
  # interpolated heatmap
  hm = Interpolations.linear_interpolation(domain, field_on_grid) # depr .LinearInterpolation(..)
  d_scalar = Vector{Float64}(undef, length(vec_preIS))

  # interpolate d_scalar for intermediate test points
  for (i, u) in enumerate(vec_preIS)
    if maximum(domain[1]) < abs(u[1]) || maximum(domain[2]) < abs(u[2])
      d_scalar[i] = 0.0
      continue
    end
    d_scalar[i] = hm(u...)
  end

  #
  weights = exp.(-d_scalar) # unscaled Gaussian
  weights ./= sum(weights)  # normalized

  # final samplable density object
  # TODO better standardize for heatmaps on manifolds
  bw = getBW(density_)[:, 1]
  @cast pts[i, j] := vec_preIS[j][i]
  bel = kde!(collect(pts), bw, weights)
  density = ManifoldKernelDensity(TranslationGroup(Ndim(bel)), bel)

  # return `<:SamplableBelief` object
  return HeatmapGridDensity(data, domain, hint_callback, bw_factor, density)
end

function Base.isapprox(
  a::HeatmapGridDensity,
  b::HeatmapGridDensity;
  atol::Real = 1e-10,
  mmd_tol::Real = 1e-2,
)
  #
  isapprox(Npts(a.densityFnc), Npts(b.densityFnc); atol) ? nothing : (return false)
  isapprox(a.densityFnc, b.densityFnc; atol = mmd_tol) ? nothing : (return false)
  isapprox(a.data, b.data; atol) ? nothing : (return false)
  isapprox(a.domain[1], b.domain[1]; atol) ? nothing : (return false)
  isapprox(a.domain[2], b.domain[2]; atol) ? nothing : (return false)

  return true
end

# legacy construct helper
function LevelSetGridNormal(
  field_on_grid::AbstractMatrix{<:Real},
  domain::Tuple{<:AbstractVector{<:Real}, <:AbstractVector{<:Real}},
  level::Real,
  sigma::Real;
  sigma_scale::Real = 3,
  hint_callback::Union{<:Function, Nothing} = nothing,
  bw_factor::Real = 0.7,  # kde spread between domain points 
  N::Int = 10000,
)
  #
  field = HeatmapGridDensity(field_on_grid, domain, hint_callback, bw_factor; N = N)
  return LevelSetGridNormal(level, sigma, float(sigma_scale), field)
end

# Field: domain (R^2/3), image (R^1/n scalar or tensor)   e.g.: x,y -> elevation ;; x, y, z, t -> EM-field (R^(4x4))
# Field( grid_x, grid_y,.... field_grid )
# Field^ = interpolator(field_at_grid, grid)
#
# FieldGrid(data_on_grid, grid_domain)  # internally does interpolation vodoo (same as Field^)
# BeliefGrid <: FieldGrid
# BeliefGrid(field_data: FieldGrid, measurement: Normal(mean: image_domain, cov: image_domain^2) ) -> domain, R_0+
# 
# calcApproxLoss(ref::BeliefGrid, appr::ManifoldKernelDensity)::Field{Real}
#  ref = Normal(ScalarField - measurement, cov)
#
