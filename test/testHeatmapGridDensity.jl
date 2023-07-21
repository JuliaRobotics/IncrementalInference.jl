
using Test
using Interpolations
using IncrementalInference
using Distributions
using TensorCast
using Distributions


##
@testset "Test HeatmapGridDensity" begin
##


x = -10:0.1:10;
y = -10:0.1:10;

img = zeros(length(x),length(y))

# lambda dispatch in this scope was REALLY slow
mv = MvNormal(Diagonal([1.0;1.0]))
g(x,y) = pdf(mv,[x;y])

for (i,x_) in enumerate(x), (j,y_) in enumerate(y)
  img[i,j] = g(x_,y_)
end

##
println("build a HeatmapGridDensity")
hgd = IIF.HeatmapGridDensity(img, (x,y), nothing, 0.07; N=1000)

@show hgd

println("test packing converters")
# check conversions to packed types
phgd = convert(PackedSamplableBelief, hgd)
hgd_ = convert(SamplableBelief, phgd)

@test isapprox( hgd, hgd_ )

# use in graph

## Check that sampling of the HMD is correct

pts_ = sample(hgd,1000)[1]

@cast pts_x[i] := pts_[i][1]
@cast pts_y[i] := pts_[i][2]

f = fit(MvNormal, hcat(pts_x,pts_y)')

@test isapprox([0;0], f.μ; atol=0.15)
@test isapprox([1 0; 0 1], f.Σ.mat; atol=0.4)

##
end