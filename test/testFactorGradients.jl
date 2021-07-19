
# PoC on jacobians for a factor

using IncrementalInference
using TensorCast
using Manifolds
using Test

##

@testset "test manual call to gradient lambda utilities" begin
##

pp = LinearRelative(MvNormal([10;0],[1 0; 0 1]))
measurement = ([10.0;0.0],)
varTypes = (ContinuousEuclid{2}, ContinuousEuclid{2})
pts = ([0;0.0], [9.5;0])

##

J__, λ_fncs, λ_sizes = IIF._prepFactorGradientLambdas(pp, measurement, varTypes, pts; h=1e-4);

##

λ_fncs[1][1]()
λ_fncs[1][2]()
λ_fncs[2][1]()
λ_fncs[2][2]()


##

J__

##

gradFct = FactorGradientsCached!(pp, measurement, varTypes, pts);

##

J = gradFct(measurement..., pts...)


##

@test norm( J - [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0] ) < 1e-4


##
end

#