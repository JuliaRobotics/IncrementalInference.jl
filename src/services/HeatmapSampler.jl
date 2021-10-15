# heatmap sampler (experimental)

@info "IncrementalInference.jl is loading tools related to Interpolations.jl."

using .Interpolations 

# only export on Requires.jl
export HeatmapDensityRegular
export sampleHeatmap, sampleLevelSetGaussian!

##


"""
    $SIGNATURES

Get the grid positions at the specified height (within the provided spreads)

DevNotes
- Recheck again that `thres` works right, as a way to discard beyond 3*sigma 
"""
function sampleHeatmap( roi::AbstractMatrix{<:Real},
                        x_grid::AbstractVector{<:Real}, 
                        y_grid::AbstractVector{<:Real};
                        sigma_scale::Real=3,
                        thres = (sigma_scale^2)  )
  #

  # truncate at sigma_scale*sigma
  _roi = thres .- roi
  mask = 0 .<= _roi

  idx2d = findall(mask)  # 2D indices
  pos = (v->[x_grid[v[1]],y_grid[v[2]]]).(idx2d)
  weights = (v->_roi[v[1],v[2]]).(idx2d)
  weights ./= sum(weights)

  # recast to the appropriate shape
  @cast kp[i,j] := pos[j][i]
  kp, weights
end


"""
    $SIGNATURES

Sample points from a regular grid scalar field level set.

Notes
- modifies first argument roi to save on memory allocations
- user should duplicate if with a deepcopy if needed
"""
function sampleLevelSetGaussian!( roi::AbstractMatrix{<:Real},
                                  sigma::Real,
                                  x_grid::AbstractVector{<:Real}, 
                                  y_grid::AbstractVector{<:Real};
                                  sigma_scale::Real=3  )
  #
  # make Gaussian
  roi .^= 2
  roi .*= 0.5/(sigma^2)
  sampleHeatmap(roi, x_grid, y_grid; sigma_scale=sigma_scale)
end


# TODO make n-dimensional, and later on-manifold
# TODO better standardize for heatmaps on manifolds
function fitKDE(support,
                weights,
                x_grid::AbstractVector{<:Real}, 
                y_grid::AbstractVector{<:Real};
                bw_factor::Real=0.7  )
  #
  # 1. set the bandwidth 
  x_spacing = Statistics.mean(diff(x_grid))
  y_spacing = Statistics.mean(diff(y_grid))
  kernel_ = bw_factor*0.5*(x_spacing + y_spacing) # 70% of the average spacing
  kernel_bw = [kernel_; kernel_]                  # same bw in x and y
  # fit KDE
  kde!(support, kernel_bw, weights)
end




function HeatmapDensityRegular( data::AbstractMatrix{<:Real}, 
                                domain::Tuple{<:AbstractVector{<:Real},<:AbstractVector{<:Real}},
                                level::Real,
                                sigma::Real;
                                sigma_scale::Real=3,
                                hist_callback::Union{<:Function, Nothing}=nothing,
                                bw_factor::Real=0.7,  # kde spread between domain points 
                                N::Int=10000  )
  #

  # select the support from raw data
  roi = data.-level
  support_, weights_ = sampleLevelSetGaussian!(roi, sigma, domain...; sigma_scale=sigma_scale)
  
  # constuct a pre-density from which to draw intermediate samples
  density_ = fitKDE(support_, weights_, domain...; bw_factor=bw_factor)
  pts_preIS, = sample(density_, N)
  
  @cast vec_preIS[j][i] := pts_preIS[i,j]
  
  # weight the intermediate samples according to interpolation of raw data
  hm = Interpolations.LinearInterpolation( domain, roi ) # interpolated heatmap
  d_scalar = Vector{Float64}( undef, length(vec_preIS) )
  
  # interpolate d_scalar for intermediate test points
  for (i,u) in enumerate(vec_preIS)
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
  bw = getBW(density_)[:,1]
  @cast pts[i,j] := vec_preIS[j][i]
  bel = kde!(collect(pts), bw, weights)
  density = ManifoldKernelDensity(TranslationGroup(Ndim(bel)), bel)

  # return `<:SamplableBelief` object
  HeatmapDensityRegular(data, domain, hist_callback, level, sigma, float(sigma_scale), bw_factor, density)
end




#