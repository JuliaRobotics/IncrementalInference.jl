# sample from weights


struct AliasingScalarSampler
  domain::Vector{Float64}
  weights::StatsBase.ProbabilityWeights
  AliasingScalarSampler() = new()
  AliasingScalarSampler(x::Vector{<:Real}, p_x::Vector{<:Real}; SNRfloor::Float64=0.0) = begin
    pxf = Float64.(p_x)
    pxf .-= quantile(pxf, SNRfloor)
    pxf[pxf .< 0.0] = 0.0
    pxf ./= norm(pxf)
    wim = StatsBase.ProbabilityWeights(pxf)
    new(x, wim)
  end
end


function rand!(ass::AliasingScalarSampler, smpls::Array{Float64})
    StatsBase.alias_sample!(ass.domain, ass.weights, smpls)
    nothing
end

function rand(ass::AliasingScalarSampler, N::Int=1)
  smpls = zeros(N)
  rand!(ass,smpls)
  smpls
end
