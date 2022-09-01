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


@testset "Packing of Rayleigh" begin
##

r = Rayleigh(1.1)
r_ = packDistribution(r)

@test r_ isa PackedSamplableBelief
@test r_ isa PackedRayleigh

r__ = unpackDistribution(r_)

@test r__ isa Rayleigh
@test isapprox(r.σ, r__.σ)

##
end


## Legacy tests

# @testset "hard-coded test of PackedPrior to Prior" begin
# ##

# pt = PackedPrior("KDE:100:[1.5015]:[-98.8276 -101.803 -98.1296 -96.1897 -99.3076 -99.1881 -101.721 -101.332 -100.431 -97.7293 -96.7652 -99.3806 -95.5593 -104.237 -94.9318 -101.691 -102.255 -98.9559 -99.3386 -99.2361 -102.483 -102.896 -97.0244 -98.9643 -99.4457 -101.298 -103.263 -2.75251 5.14065 0.327863 3.60042 -0.604114 -0.0564047 -0.804898 3.05469 1.4974 1.34657 2.22745 4.78117 1.89485 -0.48091 6.63068 0.63013 -3.11422 1.73705 5.22904 -1.73223 2.47453 1.10584 -0.0179944 3.65585 4.50016 -1.95042 98.2664 98.9983 103.748 100.789 98.4127 101.397 104.364 102.125 96.3685 103.59 99.0581 100.913 101.461 105.211 103.513 99.3325 101.201 98.05 103.508 99.9785 104.624 100.202 100.258 101.579 96.6931 95.4181 299.02 296.804 301.322 298.127 299.578 298.36 296.339 300.156 299.641 297.731 299.822 296.941 295.857 299.482 302.531 301.875 302.192 301.999 300.634 294.084 300.44]")

# tt = convert(BallTreeDensity, pt.Z)
# @test isa(tt, BallTreeDensity)

# F = Prior

# upt = convert(F, pt)

# @test isa(upt, Prior)

# ## TODO add more tests
# end



#