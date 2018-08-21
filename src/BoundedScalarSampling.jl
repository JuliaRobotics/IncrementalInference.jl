# sample from weights


struct BoundedScalarSampler
  domain::Vector{Float64}
  weights::StatsBase.ProbabilityWeights
  BoundedScalarSampler() = new()
  BoundedScalarSampler(x::Vector{<:Real}, p_x::Vector{<:Real}; SNRfloor::Float64=0.0) = begin
    pxf = Float64.(p_x)
    pxf .-= quantile(pxf, SNRfloor)
    pxf[pxf .< 0.0] = 0.0
    pxf ./= norm(pxf)
    wim = StatsBase.ProbabilityWeights(pxf)
    new(x, wim)
  end
end


function rand!(bss::BoundedScalarSampler, smpls::Array{Float64})
    StatsBase.alias_sample!(bss.domain, bss.weights, smpls)
    nothing
end

function rand(bss::BoundedScalarSampler, N::Int=1)
  smpls = zeros(N)
  rand!(bss,smpls)
  smpls
end
