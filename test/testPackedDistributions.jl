# test serialization of distributions

using IncrementalInference
using Test


##

@testset "Packing Normal" begin
##

nor = Normal()

packed = packDistribution(nor)

@test packed isa PackedSamplableBelief
@test packed isa PackedNormal

upck = unpackDistribution(packed)

@test isapprox(nor, upck)

##
end


@testset "Packing for MvNormal" begin
##


mv = MvNormal([0.0;1.0], [1 0.1; 0.1 1])
# FullNormal(
#   dim: 2
#   μ: [0.0, 1.0]
#   Σ: [1.0 0.1; 0.1 1.0]
# )




##
end


#