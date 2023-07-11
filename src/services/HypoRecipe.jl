
function Base.isapprox(
  a::HypoRecipe,
  b::HypoRecipe;
  iakws...
)
  if !(
    isapprox(a.certainidx, b.certainidx; iakws...) &&
    isapprox(a.mhidx, b.mhidx; iakws...)
  )
    @debug "HypoRecipe a vs b not the same on either .certainidx or .mhidx"
    return false
  end
  if length(a.allelements) != length(b.allelements)
    @debug "HypoRecipe different lengths on a vs b .allelements"
    return false
  end
  for (i,el) in enumerate(a.allelements)
    if !isapprox(el, b.allelements[i]; iakws...)
      @debug "HypoRecipe a vs b different on .allelements"
      return false
    end
  end
  if length(a.allelements) != length(b.allelements)
    @debug "HypoRecipe different lengths on a vs b .activehypo"
    return false
  end
  for (i,el) in enumerate(a.activehypo)
    if el[1] != b.activehypo[i][1] || !isapprox(el[2], b.activehypo[i][2]; iakws...)
      @debug "HypoRecipe a vs b different on .activehypo"
      return false
    end
  end
  return true
end

Base.:(==)(
  a::HypoRecipe, 
  b::HypoRecipe
) = isapprox(a,b)

function Base.isapprox(
  a::HypoRecipeCompute, 
  b::HypoRecipeCompute;
  iakws...
)
  if !(isnothing(a.hypotheses) && isnothing(b.hypotheses))
    return isapprox(a.hypotheses.p, b.hypotheses.p; iakws...)
  end
  if !(isnothing(a.certainhypo) && isnothing(b.certainhypo))
    return isapprox(a.certainhypo, b.certainhypo; iakws...)
  end
  
  if 0 < length(a.activehypo)
    if length(a.activehypo) == length(b.activehypo)
      return isapprox(a.activehypo, b.activehypo; iakws...)
    else
      return false
    end
  end

  return true
end


Base.:(==)(
  a::HypoRecipeCompute, 
  b::HypoRecipeCompute
) = isapprox(a,b)
