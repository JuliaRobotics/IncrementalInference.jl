# test FluxModelsDistribution and serialization

using Test
using BSON, Flux, IncrementalInference
using TensorCast

##

@testset "FluxModelsDistribution serialization" begin
##

# can a model be serialized
mdls = [Chain(Dense(5,2),Dense(2,3)); Chain(Dense(5,4), Dense(4,3))]
fxd = FluxModelsDistribution((5,),(3,),mdls,rand(5), false, false)

# check sampler is working
measd = rand(fxd, 2)
@test size( measd ) == (3,2)

# convert to packed type
fxp = convert(PackedBelief, fxd) # TODO, PackedBelief
@test fxp isa IIF.PackedFluxModelsDistribution

# convert back to hydrated object
fxu = convert(SamplableBelief, fxp)
@test fxu isa FluxModelsDistribution

measu = rand(fxu, 2)

@test measd[1] - measu[1] |> norm < 1e-6

##
end


@testset "FluxModelsDistribution serialization" begin
##

mdls = [Chain(Dense(5,2),Dense(2,3)); Chain(Dense(5,4), Dense(4,3))]
fxd = FluxModelsDistribution((5,),(3,),mdls,rand(5), false, false)

fg = initfg()

addVariable!(fg, :x0, ContinuousEuclid{3})

pr = Prior(fxd)

addFactor!(fg, [:x0;], pr)

##

smpls = sampleFactor(fg, :x0f1, 10)
@test eltype(smpls) <: AbstractVector{<:Real}
@test smpls isa Vector #{Vector{Float64}}
@test length( smpls ) == 10

##

# check local product
localProduct(fg, :x0)

solveTree!(fg)

##

saveDFG("/tmp/fg_test_flux", fg)

fg_ = loadDFG("/tmp/fg_test_flux")

ff1 = getObservation(fg_, :x0f1)
ff1.Z.shuffle[] = true

solveTree!(fg_);

# remove the testing file
Base.rm("/tmp/fg_test_flux.tar.gz")

##

end


@testset "FluxModelsDistribution as Mixture with relative factor" begin

##

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
pts_ = getBelief(fg, :x0) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test 80 < sum(-3 .< (pts) .< 3)
# at least some points should land according to the naive model
pts_ = getBelief(fg, :x1) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test 5 < sum(5 .< (pts) .< 15)


# will predict from existing fg
f1 = getObservation(fg, :x0x1f1)
predictions = map(f->f(f1.Z.components.nn.data), f1.Z.components.nn.models)


# unpack into new fg_
fg_ = loadDFG("/tmp/fg_test_flux")

# same predictions with deserialized object
f1_ = getObservation(fg_, :x0x1f1)
predictions_ = map(f->f(f1_.Z.components.nn.data), f1_.Z.components.nn.models)

# check that all predictions line up
@show norm(predictions - predictions_)
@test norm(predictions - predictions_) < 1e-6

f1_.Z.components.nn.shuffle[] = true

# test solving of the new object
solveTree!(fg_);


# prior should pin x0 pretty well
pts_ = getBelief(fg_, :x0) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test 80 < sum(-3 .< (pts) .< 3)
# at least some points should land according to the naive model
pts_ = getBelief(fg_, :x1) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test 5 < sum(5 .< (pts) .< 15)


# remove the testing file
Base.rm("/tmp/fg_test_flux.tar.gz")

##

end


@testset "MixtureFluxModels testing" begin

##

# some made up data
data = randn(10)
# Flux models
models = [Flux.Chain(softmax, Dense(10,5,σ), Dense(5,1, tanh)) for i in 1:20]
# mixture with user defined names (optional) -- could also just pass Vector or Tuple of components
mix = MixtureFluxModels(PriorCircular, models, (10,), data, (1,), 
                        (naiveNorm=Normal(),naiveUnif=Uniform()),
                        [0.7; 0.2; 0.1],
                        shuffle=false )
#

# test by add to simple graph
fg = initfg()
addVariable!(fg, :testmix, Circular)
addFactor!(fg, [:testmix;], mix)

pts = approxConv(fg, :testmixf1, :testmix);

# look at proposal distribution from the only factor on :testmix
_,pts,__, = localProduct(fg, :testmix);

saveDFG("/tmp/fg_mfx", fg)

# 
fg_ = loadDFG("/tmp/fg_mfx")

Base.rm("/tmp/fg_mfx.tar.gz")

solveTree!(fg_);

@test 10 < (getBelief(fg_, :testmix) |> getPoints |> length)

##

end





#