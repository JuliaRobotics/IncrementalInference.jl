
## Distributions to JSON/Packed types

DFG.pack(dtr::Categorical) = PackedCategorical(; p = dtr.p)
DFG.pack(dtr::Uniform) = PackedUniform(; a = dtr.a, b = dtr.b)
DFG.pack(dtr::Normal) = PackedNormal(; mu = dtr.μ, sigma = dtr.σ)
DFG.pack(dtr::ZeroMeanDiagNormal) = PackedZeroMeanDiagNormal(; diag = dtr.Σ.diag)
DFG.pack(dtr::ZeroMeanFullNormal) = PackedZeroMeanFullNormal(; cov = dtr.Σ.mat[:])
DFG.pack(dtr::DiagNormal) = PackedDiagNormal(; mu = dtr.μ, diag = dtr.Σ.diag)
DFG.pack(dtr::FullNormal) = PackedFullNormal(; mu = dtr.μ, cov = dtr.Σ.mat[:])
DFG.pack(dtr::Rayleigh) = PackedRayleigh(; sigma = dtr.σ)

## Unpack JSON/Packed to Distribution types

DFG.unpack(dtr::PackedCategorical) = Categorical(dtr.p ./ sum(dtr.p))
DFG.unpack(dtr::PackedUniform) = Uniform(dtr.a, dtr.b)
DFG.unpack(dtr::PackedNormal) = Normal(dtr.mu, dtr.sigma)
function DFG.unpack(dtr::PackedZeroMeanDiagNormal)
  return MvNormal(LinearAlgebra.Diagonal(map(abs2, sqrt.(dtr.diag))))
end # sqrt.(dtr.diag)
function DFG.unpack(dtr::PackedZeroMeanFullNormal)
  d = round(Int, sqrt(size(dtr.cov)[1]))
  return MvNormal(reshape(dtr.cov, d, d))
end
DFG.unpack(dtr::PackedDiagNormal) = MvNormal(dtr.mu, sqrt.(dtr.diag))
function DFG.unpack(dtr::PackedFullNormal)
  return MvNormal(dtr.mu, reshape(dtr.cov, length(dtr.mu), :))
end
DFG.unpack(dtr::PackedRayleigh) = Rayleigh(dtr.sigma)

