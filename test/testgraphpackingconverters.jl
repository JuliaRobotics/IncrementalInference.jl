# packing and unpacking of Graph related types

using IncrementalInference
using DistributedFactorGraphs

using Test

dfg = GraphsDFG{SolverParams}(params=SolverParams())
setSerializationModule!(dfg, Main)

@testset "hard-coded test of PackedPrior to Prior" begin

pt = PackedPrior("KDE:100:[1.5015]:[-98.8276 -101.803 -98.1296 -96.1897 -99.3076 -99.1881 -101.721 -101.332 -100.431 -97.7293 -96.7652 -99.3806 -95.5593 -104.237 -94.9318 -101.691 -102.255 -98.9559 -99.3386 -99.2361 -102.483 -102.896 -97.0244 -98.9643 -99.4457 -101.298 -103.263 -2.75251 5.14065 0.327863 3.60042 -0.604114 -0.0564047 -0.804898 3.05469 1.4974 1.34657 2.22745 4.78117 1.89485 -0.48091 6.63068 0.63013 -3.11422 1.73705 5.22904 -1.73223 2.47453 1.10584 -0.0179944 3.65585 4.50016 -1.95042 98.2664 98.9983 103.748 100.789 98.4127 101.397 104.364 102.125 96.3685 103.59 99.0581 100.913 101.461 105.211 103.513 99.3325 101.201 98.05 103.508 99.9785 104.624 100.202 100.258 101.579 96.6931 95.4181 299.02 296.804 301.322 298.127 299.578 298.36 296.339 300.156 299.641 297.731 299.822 296.941 295.857 299.482 302.531 301.875 302.192 301.999 300.634 294.084 300.44]")

tt = convert(BallTreeDensity, pt.Z)
@test isa(tt, BallTreeDensity)

F = Prior

upt = convert(F, pt)

@test isa(upt, Prior)

  # TODO add more tests
end



fg = initfg()
N=100

doors = zeros(1,4)
doors[1,1:4] = [-100.0;0.0;100.0;300.0]
pd = kde!(doors,[3.0])
pd = resample(pd,N);
bws = getBW(pd)[:,1]
bws2 = Array{Float64,2}(undef, length(bws),1)
bws2[:,1] = bws[:]
doors2 = getPoints(pd);
v1 = addVariable!(fg,:x1, ContinuousScalar,N=N)
f1  = addFactor!(fg,[:x1], Prior(kde!(doors2, bws2[:])))

v2 = addVariable!(fg,:x2, ContinuousScalar, N=N)
lc = LinearConditional(  Normal(50.0, 2.0) )
f2 = addFactor!(fg, [:x1; :x2], lc)


@testset "Testing conversion to packed function node data structure and back" begin

topack = solverData(f1)
dd = convert(PackedFunctionNodeData{PackedPrior},topack)
upd = convert(FunctionNodeData{CommonConvWrapper{Prior}}, dd)

@test compare(topack, upd)

topack = solverData(f2)
dd =  convert(IncrementalInference.PackedFunctionNodeData{PackedLinearConditional},topack)
upd = convert(IncrementalInference.FunctionNodeData{CommonConvWrapper{LinearConditional}}, dd)

@test compare(topack, upd)

end

@testset "Testing conversion to packed variable node data structure and back" begin

dat = solverData(getVariable(fg,:x1))

# dat.BayesNetVertID

pd = pack(dfg, dat)
unpckd = unpack(dfg, pd)

@test compareFields(dat, unpckd, skip=[:softtype])
@test compareFields(dat.softtype, unpckd.softtype)
@test isa(dat.softtype, ContinuousScalar)
@test isa(unpckd.softtype, ContinuousScalar)

end





#
