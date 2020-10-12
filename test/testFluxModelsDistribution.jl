# test FluxModelsDistribution and serialization

using Test
using BSON, Flux, IncrementalInference



@testset "FluxModelsDistribution serialization" begin

# can a model be serialized
mdls = [Chain(Dense(5,2),Dense(2,3)); Chain(Dense(5,4), Dense(4,3))]
fxd = FluxModelsDistribution((5,),(3,),mdls,rand(5), false, false)

# check sampler is working
measd = rand(fxd, 2)
@test measd |> size == (3,2)

# convert to flat string
fxp = convert(PackedSamplableBelief, fxd)
@test fxp isa String

# convert back to hydrated object
fxu = convert(SamplableBelief, fxp)
@test fxu isa FluxModelsDistribution

measu = rand(fxu, 2)

@test measd[1] - measu[1] |> norm < 1e-6


end




@testset "FluxModelsDistribution serialization" begin



mdls = [Chain(Dense(5,2),Dense(2,3)); Chain(Dense(5,4), Dense(4,3))]
fxd = FluxModelsDistribution((5,),(3,),mdls,rand(5), false, false)

fg = initfg()

addVariable!(fg, :x0, ContinuousEuclid{3})

pr = Prior(fxd)

addFactor!(fg, [:x0;], pr)

# check local product
localProduct(fg, :x0)

solveTree!(fg)


saveDFG("/tmp/fg_test_flux", fg)

fg_ = loadDFG("/tmp/fg_test_flux")

ff1 = getFactorType(fg_, :x0f1)
ff1.Z.shuffle[] = true

solveTree!(fg_);



end

#