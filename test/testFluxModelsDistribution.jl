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

# remove the testing file
Base.rm("/tmp/fg_test_flux.tar.gz")

end




@testset "FluxModelsDistribution as Mixture with relative factor" begin


mdls = [Chain(Dense(10,50, relu),Dense(50,20),softmax, Dense(20,1, tanh)) for i in 1:50];
fxd = FluxModelsDistribution((10,),(1,),mdls,rand(10), false, false)


fg = initfg()

addVariable!(fg, :x0, ContinuousEuclid{1})
addVariable!(fg, :x1, ContinuousEuclid{1})

# a prior
pr = Prior(Normal())
addFactor!(fg, [:x0;], pr)

# a relative mixture network
mfx = Mixture(LinearRelative, (naive=Normal(10, 10), nn=fxd), [0.5;0.5])
addFactor!(fg, [:x0;:x1], mfx)

# and test overall serialization before solving
saveDFG("/tmp/fg_test_flux", fg)

# solve existing fg
solveTree!(fg)

# prior should pin x0 pretty well
@test 80 < sum(-3 .< (getBelief(fg, :x0) |> getPoints) .< 3)
# at least some points should land according to the naive model
@test 5 < sum(5 .< (getBelief(fg, :x1) |> getPoints) .< 15)


# will predict from existing fg
f1 = getFactorType(fg, :x0x1f1)
predictions = map(f->f(f1.components.nn.data), f1.components.nn.models)


# unpack into new fg_
fg_ = loadDFG("/tmp/fg_test_flux")

# same predictions with deserialized object
f1_ = getFactorType(fg_, :x0x1f1)
predictions_ = map(f->f(f1_.components.nn.data), f1_.components.nn.models)

# check that all predictions line up
@show norm(predictions - predictions_)
@test norm(predictions - predictions_) < 1e-6

f1_.components.nn.shuffle[] = true

# test solving of the new object
solveTree!(fg_);


# prior should pin x0 pretty well
@test 80 < sum(-3 .< (getBelief(fg_, :x0) |> getPoints) .< 3)
# at least some points should land according to the naive model
@test 5 < sum(5 .< (getBelief(fg_, :x1) |> getPoints) .< 15)


# remove the testing file
Base.rm("/tmp/fg_test_flux.tar.gz")

end




#