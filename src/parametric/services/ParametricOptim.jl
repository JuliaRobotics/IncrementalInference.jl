
# experimental
function optimizeManifold_FD(
  M::AbstractManifold, 
  cost::Function,
  x0::AbstractArray;
  algorithm = Optim.ConjugateGradient(; manifold=ManifoldWrapper(M))
)
  # finitediff setup
  r_backend = ManifoldDiff.TangentDiffBackend(
    ManifoldDiff.FiniteDifferencesBackend()
  )
  
  ## finitediff gradient (non-manual)
  function costgrad_FD!(X,p)
    X .= ManifoldDiff.gradient(M, cost, p, r_backend)
    X
  end

  Optim.optimize(cost, costgrad_FD!, x0, algorithm)
end
