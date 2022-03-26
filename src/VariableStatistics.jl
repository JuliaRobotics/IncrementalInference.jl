


Statistics.mean(vartype::InferenceVariable, args...; kwargs...) = mean(getManifold(vartype), args...; kwargs...)
Statistics.std(vartype::InferenceVariable, args...; kwargs...) = std(getManifold(vartype), args...; kwargs...)
Statistics.var(vartype::InferenceVariable, args...; kwargs...) = var(getManifold(vartype), args...; kwargs...)

function Statistics.cov(vartype::InferenceVariable, ptsArr::AbstractVector; basis::Manifolds.AbstractBasis = Manifolds.DefaultOrthogonalBasis(), kwargs...)
  return cov(getManifold(vartype), ptsArr; basis, kwargs... )
end

# To replace calcCovarianceBasic
function calcStdBasicSpread(vartype::InferenceVariable, ptsArr::Vector{P}) where P
  σ = std(vartype, ptsArr)

  #if no std yet, set to 1
  msst = 1e-10 < σ ? σ : 1.0
  return msst
end

#TODO consolidate
function calcMeanCovar(vari::DFGVariable, solvekey=:default)

  pts = getSolverData(vari, solvekey).val
  μ = mean(getManifold(vari), pts)
  Σ = cov(getVariableType(vari), pts)
  return μ, Σ
end


