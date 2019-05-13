import Base.convert

"""
Converter: Packed_MvNormal -> MvNormal
"""
function convert(::Type{Distributions.MvNormal}, pv::Dict{String, Any})
    len = length(Float64.(pv["mean"]))
    mat = reshape(Float64.(pv["cov"]), len, len)
    return Distributions.MvNormal(Float64.(pv["mean"]), mat)
end

"""
Converter: MvNormal -> Packed_MvNormal
"""
function convert(::Type{Dict{String, Any}}, mvNormal::Distributions.MvNormal)
    v = mvNormal.Σ.mat[:]
    # TODO: Confirm that the Dict{String, Any} is necessary, could just go string and back
    return JSON2.read(JSON2.write(Packed_MvNormal(mvNormal.μ, v, "MvNormal")), Dict{String,Any})
end

"""
Converter: Packed_Normal -> Normal
"""
function convert(::Type{Distributions.Normal}, pv::Dict{String, Any})
    return Distributions.Normal(Float64(pv["mean"]), Float64(pv["std"]))
end

"""
Converter: Normal -> Packed_Normal
"""
function convert(::Type{Dict{String, Any}}, normal::Distributions.Normal)
    return JSON2.read(JSON2.write(Packed_Normal(normal.μ, normal.σ, "Normal")), Dict{String,Any})
end


"""
Converter: Packed_AliasingScalarSampler -> AliasingScalarSampler
"""
function convert(::Type{AliasingScalarSampler}, pv::Dict{String, Any})
    snrFloor = haskey(pv, "quantile") && pv["quantile"] != nothing ? Float64(pv["quantile"]) : 0.0
    sampler = AliasingScalarSampler(Float64.(pv["samples"]), Float64.(pv["weights"]); SNRfloor = snrFloor)
    return sampler
end

"""
Converter: AliasingScalarSampler -> Packed_AliasingScalarSampler
"""
function convert(::Type{Dict{String, Any}}, sampler::AliasingScalarSampler)
    packed = Packed_AliasingScalarSampler(sampler.domain, sampler.weights.values, 0.0, "AliasingScalarSampler")
    return JSON2.read(JSON2.write(packed), Dict{String,Any})
end

# Catch-all with error
function convert(::Type{Dict{String, Any}}, distribution::T) where T <: Sampleable
    error("Converter does not exist for type $(T)->Dict{String, Any}. Please introduce one in IncrementalInference.")
end

function convert(::Type{T}, packed::Dict{String, Any}) where T <: Sampleable
    error("Converter does not exist for type Dict{String, Any}->$(T). Please introduce one in IncrementalInference.")
end
