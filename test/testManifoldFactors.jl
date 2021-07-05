

using IncrementalInference
using Manifolds
using LinearAlgebra
using StaticArrays

##

@testset "test manifold factors" begin

##

f = ManifoldFactor(SpecialOrthogonal(3), MvNormal([0.1, 0.02, 0.01]))
s = getSample(f,10)[1]
s[1]

f = ManifoldFactor(SpecialEuclidean(2), MvNormal([0.1, 0.2, 0.01]))
s = getSample(f,10)[1]
s[1]


f = ManifoldPrior(SpecialOrthogonal(2), MvNormal([0.1]), SA[1.0 0; 0 1])
meas = getSample(f,10)[1]
meas[1]
f.(meas, Ref(SA[1.0 0; 0 1]))

f = ManifoldPrior(SpecialOrthogonal(3), MvNormal([0.1, 0.02, 0.01]), SA[1.0 0 0; 0 1 0; 0 0 1])
s = getSample(f,10)[1]
s[1]

f = ManifoldPrior(SpecialEuclidean(2), MvNormal([0.1, 0.2, 0.01]), ProductRepr(SA[0,0], SA[1.0 0; 0 1]))
s = getSample(f,10)[1]
s[1]

##

end

##