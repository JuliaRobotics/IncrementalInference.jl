# test nullhypo with n-dim partial

using Test
using IncrementalInference

##


@testset "test n-dimensional partial" begin

##

fg = initfg()
addVariable!(fg, :x0, ContinuousEuclid{3})

addFactor!(fg, [:x0;], PartialPrior(ContinuousEuclid{3}, MvNormal(zeros(2), ones(2)), (2,3)) )

addVariable!(fg, :x1, ContinuousEuclid{3})
addFactor!(fg, [:x1;], PartialPrior(ContinuousEuclid{3}, Normal(10,1),(1,)))
addFactor!(fg, [:x0; :x1], LinearRelative(MvNormal([10;0;0.0], ones(3))) )

##

initAll!(fg)

##

destlbl = :x0

dens = Vector{ManifoldKernelDensity}()
factors = getFactor.(fg, ls(fg, destlbl))
inferdim = IIF.proposalbeliefs!(fg, destlbl, factors, dens )

oldBel = getBelief(fg, destlbl)
oldpts = getPoints(oldBel)

varType = getVariableType(fg, destlbl)
pGM = getPoints( AMP.manifoldProduct(dens, getManifold(varType), N=100, oldPoints=oldpts), false )
# pGM = AMP.productbelief(oldpts, getManifold(varType), dens, 100, asPartial=false )


##

densPts, inferdim = propagateBelief(fg, :x0, :, needFreshMeasurements=true )

##

solveTree!(fg);

##

@warn "WIP on testPartialNH.jl during transition to Manifolds.jl"
@test isapprox( calcMeanMaxSuggested(fg, :x0, :default).suggested, [0;0;0], atol=1)
@test isapprox( calcMeanMaxSuggested(fg, :x1, :default).suggested, [10;0;0], atol=1)

##

end


@testset "test n-dimensional partial with nullhypo" begin

##

fg = initfg()
addVariable!(fg, :x0, ContinuousEuclid{3})

addFactor!(fg, [:x0;], PartialPrior(ContinuousEuclid{3}, MvNormal(zeros(2), ones(2)),(2,3)) , nullhypo=0.2)

addVariable!(fg, :x1, ContinuousEuclid{3})
addFactor!(fg, [:x1;], PartialPrior(ContinuousEuclid{3}, Normal(10,1),(1,)))
addFactor!(fg, [:x0; :x1], LinearRelative(MvNormal([10;0;0.0], ones(3))) , nullhypo=0.2)

##

solveTree!(fg);

##

@warn "WIP testPartialNH.jl during transition to Manifolds.jl"
@test isapprox( calcMeanMaxSuggested(fg, :x0, :default).suggested, [0;0;0], atol=1)
@test isapprox( calcMeanMaxSuggested(fg, :x1, :default).suggested, [10;0;0], atol=2)

##

end

#