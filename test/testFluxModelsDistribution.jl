# test FluxModelsDistribution and serialization

using Test
using BSON, Flux, IncrementalInference



@testset "FluxModelsDistribution serialization" begin

# can a model be serialized
mdls = [Chain(Dense(5,2),Dense(2,3)); Chain(Dense(5,4), Dense(4,3))]
fxd = FluxModelsDistribution((5,),(3,),mdls,rand(5), false, false)

# convert to flat string
fxp = convert(PackedSamplableBelief, fxd)
@test fxp isa String

# convert back to hydrated object
fxu = convert(SamplableBelief, fxp)
@test fxu isa FluxModelsDistribution

measd = IIF.rand(fxd, 2)
measu = IIF.rand(fxu, 2)

@test measd[1] - measu[1] |> norm < 1e-6


end




@testset "FluxModelsDistribution serialization" begin


mdls = [Chain(Dense(5,2),Dense(2,3)); Chain(Dense(5,4), Dense(4,3))]
fxd = FluxModelsDistribution((5,),(3,),mdls,rand(5), false, false)



fg = initfg()

addVariable!(fg, :x0, ContinuousEuclid{3})

pr = Prior(fxd)

addFactor!(fg, [:x0;], pr)


end

#