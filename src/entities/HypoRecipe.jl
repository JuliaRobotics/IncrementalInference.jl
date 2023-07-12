


Base.@kwdef struct HypoRecipe
  certainidx::Vector{Int}
  allelements::Vector{Vector{Int}}
  activehypo::Vector{Tuple{Int,Vector{Int}}}
  mhidx::Vector{Int}
end

Base.@kwdef struct HypoRecipeCompute{
  HP <: Union{Nothing, <:Distributions.Categorical{Float64, <:AbstractVector{Float64}}},
  CH <: Union{Nothing, <:AbstractVector{<:Integer}},
}
  """ multi hypothesis settings #NOTE no need for a parameter as type is known from `parseusermultihypo` """
  hypotheses::HP = nothing
  """ categorical to select which hypothesis is being considered during convolution operation """
  certainhypo::CH = nothing
  """ subsection indices to select which params should be used for this hypothesis evaluation """
  activehypo::Vector{Int} = Int[]
end

