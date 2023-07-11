


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


function Base.isapprox(
  a::HypoRecipe, 
  b::HypoRecipe
)
  if !(isnothing(a.hypotheses) && isnothing(B.hypotheses))
    return isapprox(a.hypotheses.p, b.hypotheses.p)
  end
  if !(isnothing(a.certainhypo) && isnothing(B.certainhypo))
    return isapprox(a.certainhypo, b.certainhypo)
  end
  
  if 0 < length(a.activehypo)
    if length(a.activehypo) == length(b.activehypo)
      return isapprox(a.activehypo, b.activehypo)
    else
      return false
    end
  end

  return true
end


Base.==(
  a::HypoRecipe, 
  b::HypoRecipe
) = isapprox(a,b)
