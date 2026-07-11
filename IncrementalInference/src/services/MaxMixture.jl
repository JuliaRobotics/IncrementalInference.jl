
## ================================================================================================
## Experimental specialized dispatch for Mixture
## ================================================================================================
# To sort out how to dispatch on specialized functions.
# related to #931 and #1069

struct MaxMixture <: AbstractMaxMixtureSolver
  p::Vector{Float64}
  # the chosen component to be used for the optimization
  choice::Base.RefValue{Int}
end

function getMeasurementParametric(s::Mixture)
  meas = map(c -> getMeasurementParametric(c)[1], values(s.components))
  iΣ = map(c -> getMeasurementParametric(c)[2], values(s.components))
  return meas, iΣ
end

function _calcFactorMahalanobis(cfp, meas, iΣ, variables...)
  res = cfp.calcfactor!(meas, variables...)
  r = res' * iΣ * res
  return r
end

# DEV NOTE: function with other options including select once and use
# function (cfp::CalcFactorMahalanobis{<:CalcFactor, MaxMixture})(variables...)
#   if cfp.specialAlg.choice[] == 0
#     #calculate all mixture options
#     r = [_calcFactorMahalanobis(cfp, cfp.meas[i], cfp.iΣ[i], variables...) for i = 1:length(cfp.meas)]

#     p = cfp.specialAlg.p

#     k = size(cfp.iΣ[1], 2)
#     # α = 1 ./ sqrt.(2pi .* k .* det.(inv.(cfp.iΣ)))
#     α = sqrt.(det.(cfp.iΣ) ./ ((2pi)^k))

#     # mm, at = findmax(α .* p .* exp.(-0.5 .* r))
#     # mm = sum(α .* p .* exp.(-0.5 .* r) )

#     mm, at = findmin( 0.5 .* r .- log.(α .* p))
#     # mm = -log(sum(α .* p .* exp.(-0.5 .* r) ))
#     # return mm + maximum(log.(α .* p))

#     cfp.specialAlg.choice[] = at

#     return r[at] 

#   else
#     at = cfp.specialAlg.choice[]
#     return _calcFactorMahalanobis(cfp, cfp.meas[at], cfp.iΣ[at], variables...)
#   end

# end

# function (cfp::CalcFactorMahalanobis{FT, N, C, MEAS, D, L, MaxMixture})(variables...) where {FT, N, C, MEAS, D, L}
#   r = [
#     _calcFactorMahalanobis(cfp, cfp.meas[i], cfp.iΣ[i], variables...) for
#     i = 1:length(cfp.meas)
#   ]

#   p = cfp.specialAlg.p

#   k = size(cfp.iΣ[1], 2)
#   # α = 1 ./ sqrt.(2pi .* k .* det.(inv.(cfp.iΣ)))
#   α = sqrt.(det.(cfp.iΣ) ./ ((2pi)^k))

#   mm, at = findmin(r .- log.(α .* p))
#   # mm = -log(sum(α .* p .* exp.(-0.5 .* r) ))
#   return mm + maximum(log.(α .* p))
# end

## ================================================================================================
## Experimental specialised dispatch for multihypo and nullhypo
## ================================================================================================
#TODO better dispatch

struct MaxMultihypo <: AbstractMaxMixtureSolver
  multihypo::Vector{Float64}
end
struct MaxNullhypo <: AbstractMaxMixtureSolver
  nullhypo::Float64
end

# function (cfp::CalcFactorMahalanobis{FT, N, C, MEAS, D, L, Nothing})(X1, L1, L2) where {FT, N, C, MEAS, D, L}
#   mh = cfp.specialAlg.multihypo
#   @assert length(mh) == 3 "multihypo $mh  not supported with parametric, length should be 3"
#   @assert mh[1] == 0 "multihypo $mh  not supported with parametric, first should be 0"

#   #calculate both multihypo options
#   r1 = cfp(X1, L1)
#   r2 = cfp(X1, L2)
#   r = [r1, r2]

#   # hacky multihypo to start of with 
#   mm, at = findmin(r .* (1 .- mh[2:end]))
#   nat = at == 1 ? 1 : 2
#   k = length(X1) * one(r1) * 1e-3
#   return r[at] + r[nat] * k
# end

# function (cfp::CalcFactorMahalanobis{FT, N, C, MEAS, D, L, MaxNullhypo})(X1, X2) where {FT, N, C, MEAS, D, L}
#   nh = cfp.specialAlg.nullhypo
#   @assert nh > 0 "nullhypo $nh not as expected"

#   #calculate factor residual
#   res = cfp(cfp.meas[1], X1, X2)
#   r1 = res' * cfp.iΣ * res

#   # compare to uniform nullhypo
#   r2 = length(res) * one(r1)
#   r = [r1, r2]
#   mm, at = findmin(r .* [nh, (1 - nh)])

#   residual = at == 1 ? r1 : r1 * 1e-3

#   return residual

#   # rand residual option
#   # idx = rand(Categorical([(1-nh), nh]))
#   # nh == 0.05 && cfp.varOrder==[:x1,:l1] && println("$idx -> $(r1.value), $r2")
#   # return r[idx] 

# end
