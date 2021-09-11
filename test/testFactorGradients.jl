
# PoC on jacobians for a factor

using IncrementalInference
using TensorCast
using Manifolds
using Test

# overloading with new dispatch
import IncrementalInference: getSample, getManifold

##

@testset "test manual call to gradient lambda utilities" begin
##

pp = LinearRelative(MvNormal([10;0],[1 0; 0 1]))
measurement = [10.0;0.0]
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

@test norm( J__ - [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0] ) < 1e-4


## build new functor container

gradFct = FactorGradientsCached!(pp, varTypes, measurement, pts);

## test grad calc for current values

J_c = gradFct()

@test norm( J_c - [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0] ) < 1e-4

##

J = gradFct(measurement, pts...)


##

@test norm( J - [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0] ) < 1e-4

## check on transmitted info per coords

ret = calcPerturbationFromVariable(gradFct, [1=>[1;1]])

# the fromVar itself should be zero
@test length(ret[1]) == 2
@test isapprox( ret[1], [0;0], atol=1e-6 )

# the other variable
@test length(ret[2]) == 2
@test isapprox( ret[2], [1;1], atol=1e-6 )

##
end

##

struct TestPartialRelative2D{B <: SamplableBelief} <: IIF.AbstractRelativeMinimize
  Z::B
  partial::Tuple{Int}
end
# standard helper with partial set
TestPartialRelative2D(z::SamplableBelief) = TestPartialRelative2D(z, (2,))

# imported earlier for overload
getManifold(fnc::TestPartialRelative2D) = TranslationGroup(2)
getSample(cf::CalcFactor{<:TestPartialRelative2D}) = rand(cf.factor.Z, 1)

# currently requires residual to be returned as a tangent vector element
(cf::CalcFactor{<:TestPartialRelative2D})(z, x1, x2) = x2[2:2] - (x1[2:2] + z[1:1])

##

@testset "test a partial, binary relative factor perturbation (a new user factor)" begin
##


tpr = TestPartialRelative2D(Normal(10,1))

measurement = [10.0;]
pts = ([0;0.0], [0;10.0])
varTypes = (ContinuousEuclid{2}, ContinuousEuclid{2})

## construct the lambdas
gradients = FactorGradientsCached!(tpr, varTypes, measurement, pts);

## calculate new gradients
J = gradients(measurement, pts...)


## check on transmitted info per coords

ret = calcPerturbationFromVariable(gradients, [1=>[1;1]])

# the fromVar itself should be zero
@test length(ret[1]) == 2
@test isapprox( ret[1], [0;0], atol=1e-6 )

# the other variable only affects the second coordinate dimension
@test length(ret[2]) == 2
@test isapprox( ret[2], [0;1], atol=1e-6 )

## check the reverse perturbation

ret = calcPerturbationFromVariable(gradients, [2=>[1;1]])

# only the first coordinate dimension is affected
@test length(ret[1]) == 2
@test isapprox( ret[1], [0;1], atol=1e-6 )

## the fromVar itself should be zero
@test length(ret[2]) == 2
@test isapprox( ret[2], [0;0], atol=1e-6 )

##
end


#