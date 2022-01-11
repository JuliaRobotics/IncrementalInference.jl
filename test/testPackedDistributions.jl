# test serialization of distributions

using IncrementalInference
using Test


##

@testset "Packing Categorical" begin
##

ctg = Categorical(5)

packed = packDistribution(ctg)

@test packed isa PackedSamplableBelief
@test packed isa IncrementalInference.PackedCategorical

upck = unpackDistribution(packed)

@test upck isa Categorical
@test isapprox(ctg, upck)

##
end


@testset "Packing Normal" begin
##

nor = Normal()

packed = packDistribution(nor)

@test packed isa PackedSamplableBelief
@test packed isa PackedNormal

upck = unpackDistribution(packed)

@test upck isa Normal
@test isapprox(nor, upck)

##
end


@testset "Packing for FullNormal" begin
##

mv = MvNormal([0.0;1.0], [1 0.1; 0.1 1])
packed = packDistribution(mv)

@test packed isa PackedSamplableBelief
@test packed isa PackedFullNormal

upck = unpackDistribution(packed)

@test upck isa FullNormal
@test isapprox(mv, upck)

##
end


@testset "Packing for DiagNormal" begin
##

mv = MvNormal([0.0; 1.0], [4.0; 4.0])
packed = packDistribution(mv)

@test packed isa PackedSamplableBelief
@test packed isa PackedDiagNormal

upck = unpackDistribution(packed)

@test upck isa DiagNormal
@test isapprox(mv, upck)

##
end

@testset "Packing for ZeroMeanFullNormal" begin
##

mv = MvNormal([1 0.1; 0.1 1])
packed = packDistribution(mv)

@test packed isa PackedSamplableBelief
@test packed isa PackedZeroMeanFullNormal

upck = unpackDistribution(packed)

@test upck isa ZeroMeanFullNormal
@test isapprox(mv, upck)

##
end


@testset "Packing for ZeroMeanDiagNormal" begin
##

mv = MvNormal([4.0;4.0])
packed = packDistribution(mv)

@test packed isa PackedSamplableBelief
@test packed isa PackedZeroMeanDiagNormal

upck = unpackDistribution(packed)

@test upck isa ZeroMeanDiagNormal
@test isapprox( mv.Σ.diag, upck.Σ.diag )
# @test isapprox(mv, upck)

##
end



@testset "Packing for AliasingScalarSampler" begin
##

bss = AliasingScalarSampler([1.0;2.0], [0.6;0.4])
packed = packDistribution(bss)

@test packed isa PackedSamplableBelief
@test packed isa IncrementalInference.PackedAliasingScalarSampler

upck = unpackDistribution(packed)

@test upck isa AliasingScalarSampler
@test isapprox( bss.domain, upck.domain )
@test isapprox( bss.weights.values, upck.weights.values )
@test isapprox( bss.weights.sum, upck.weights.sum )

##
end




@testset "Packing for ManifoldKernelDensity" begin
##

T = ContinuousEuclid{2}
pts = [randn(2) for _ in 1:50]

mkd = manikde!(T, pts)
packed = packDistribution(mkd)

@test packed isa PackedSamplableBelief
@test packed isa IncrementalInference.PackedManifoldKernelDensity

upck = unpackDistribution(packed)

@test upck isa ManifoldKernelDensity

@test isapprox(mkd, upck)

##
end




#