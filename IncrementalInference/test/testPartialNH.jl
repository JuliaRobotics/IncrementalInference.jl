# test nullhypo with n-dim partial

using Test
using IncrementalInference

##


@testset "test 3-dimensional partial on two variable graph w (1,)--(2,3)" begin

##

  fg = initfg()
  addVariable!(fg, :x0, ContinuousEuclid{3})

  addFactor!(fg, [:x0;], PartialPrior(ContinuousEuclid{3}, MvNormal(zeros(2), ones(2)), (2,3)) )

  addVariable!(fg, :x1, ContinuousEuclid{3})
  addFactor!(fg, [:x1;], PartialPrior(ContinuousEuclid{3}, Normal(10,1),(1,)))
  addFactor!(fg, [:x0; :x1], LinearRelative(MvNormal([10;0;0.0], ones(3))) )

## check observability field is properly computed during proposals and product
  vsym = :x0
  useinitfct = IncrementalInference.listFactors_Initialized(fg, vsym; _neighbors = lsf(fg, vsym))
  dens = Vector{HomotopyDensityLive}()
  _observability = IncrementalInference.proposalbeliefs!(fg, vsym, map(x -> getFactor(fg, x), useinitfct), dens)

  @test dens[1].observability == [0;1;1]
  @test _observability == [0;1;1]

  # take the product
  hode = manifoldProduct(
    dens;
    MC = 1,
  )

  @test hode.observability == [0;1;1]

##

doautoinit!(fg, [:x0;]; singles = true)

@test getBelief(fg, :x0).observability == [0;1;1]


##

initAll!(fg)


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