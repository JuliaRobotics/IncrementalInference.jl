# sample from weights
# import IncrementalInference: AliasingScalarSampler

"""
$(TYPEDEF)

Sampler from intensity map given Euclidean domain `x` and probability weights `p_x`.

Example

`AliasingScalarSampler(x::Vector{<:Real}, p_x::Vector{<:Real}; SNRfloor::Float64=0.0)`
"""
struct AliasingScalarSampler
  domain::Vector{Float64}
  weights::StatsBase.ProbabilityWeights
  AliasingScalarSampler() = new()
  AliasingScalarSampler(x::Vector{<:Real}, p_x::Vector{<:Real}; SNRfloor::Float64=0.0) = begin

    # pxf = Float64.(p_x)
    # pxf .-= quantile(pxf, SNRfloor)
    # pxf[pxf .< 0.0] = 0.0
    # pxf ./= norm(pxf)
    # new implementation : pxf should be an empirical pmf before use in statsbase
    pxf = Float64.(p_x)
    pxf[pxf.<0.0] .= 0.0 # no negative values!
    pxf ./=sum(pxf)  # must sum to 1
        pxf2 = pxf .- quantile(pxf,SNRfloor) # remove lowest quantile
        pxf2[pxf2.<0.0] .= 0.0
        pxf2s = sum(pxf2)
        pxf[:] = 1e-10 < pxf2s ? pxf2 : pxf
        pxf ./= sum(pxf)
        sum(isnan.(pxf)) == 0 ? nothing : error("AliasingScalarSampler got NaN because of particular values in p_x")
        # pxf .-= quantile(pxf,SNRfloor) # remove lowest quantile
        # pxf[pxf.<0.0] .= 0.0
        # pxf ./=sum(pxf)
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


function convert(::Type{AliasingScalarSampler}, str::AS) where {AS <: AbstractString}
  sstr = split(split(str, "AliasingScalarSampler")[end], '[')
  sstr = sstr[length.(sstr) .> 2]
  ssstr = split.(sstr, ']')
  domain = parse.(Float64, strip.(split(ssstr[1][1], ',')) )
  weight = parse.(Float64, strip.(split(ssstr[2][1], ',')) )
  AliasingScalarSampler(domain, weight)
end
# str = "IncrementalInference.AliasingScalarSampler([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0,24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0], [0.156102, 0.117163, 0.221049, 0.275905, 0.0488494, 0.0731541, 0.107584,0.313848, 0.0309002, 0.0, 0.0384554, 0.155308, 0.276917, 0.0271168, 0.293263, 0.171316, 0.27459, 0.175323, 0.0535772, 0.181663, 0.295042, 0.104593, 0.0472137, 0.326016, 0.055283, 0.0737767, 0.302647, 0.0291257, 0.0206642, 0.223375])"
