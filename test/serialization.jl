using Test
using IncrementalInference

# Test that we can serialize and deserialize these
const PackedType = Dict{String, Any}

@testset "Distribution Serialization" begin

# Distributions that support serialization
supportedDistributions =
    [Normal(5.0),
     MvNormal([0, 0, 0, 0], Matrix(Diagonal(ones(4)))),
     AliasingScalarSampler(ones(100), round.(rand(100), digits=3); )
     # BallTreeDensity
     ]

a = convert(Dict{String, Any}, supportedDistributions[1])
b = convert(_evalType(a["distType"]), a)
supportedDistributions[1] == b

for d in supportedDistributions
    @info " - Testing $(typeof(d))"
    j = convert(PackedType, d)
    !haskey(j, "distType") && error("$(typeof(d)) when serialized does not contain the required field 'distType' - please make sure this is included.")
    dBack = convert(_evalType(j["distType"]), j)
    jBack = convert(PackedType, dBack)
    @test compareAll(j, jBack)
    # @test d == dBack
end

# Test that when we cannot, we fail correctly d = Binomial(10, 0.6)
@test_throws Exception convert(PackedType, d)
@test_throws Exception convert(Binomial, Dict{String, Any}())

end
