# packing and unpacking of Graph related types

using IncrementalInference
using DistributedFactorGraphs
using Manifolds

using Test

##

@testset "Serialization of SamplableBelief types" begin
##

td = Uniform()

ptd = pack(td) # TODO, PackedBelief
utd = unpack(td)

@test td.a - utd.a |> abs < 1e-10
@test td.b - utd.b |> abs < 1e-10

##
end
##


fg = initfg()
N=100

# doors = zeros(1,4)
T = ContinuousScalar
doors = [[-100.0;],[0.0;],[100.0;],[300.0;]]
pd = manikde!(T, doors; bw=[3.0;])
pd = resample(pd,N);
bws = getBW(pd)[:,1]
doors2 = getPoints(pd);
v1 = addVariable!(fg,:x1, ContinuousScalar,N=N)
f1  = addFactor!(fg,[:x1], Prior(manikde!(T, doors2; bw=bws)))

v2 = addVariable!(fg,:x2, ContinuousScalar, N=N)
lc = LinearRelative(  Normal(50.0, 2.0) )
f2 = addFactor!(fg, [:x1; :x2], lc)

##

@testset "Testing conversion to packed function node data structure and back" begin
##

topack = DFG.getObservation(f1)
dd = DFG.pack(topack)
upd = DFG.unpack(dd)

@test topack == upd

topack = DFG.getObservation(f2)
dd =  DFG.pack(topack)
upd = DFG.unpack(dd)

@test topack == upd

##
end

@testset "Testing conversion to packed variable node data structure and back" begin
##

dat = getState(getVariable(fg,:x1), :default)

# dat.BayesNetVertID

pd = packState(dat)
unpckd = unpackState(pd)

@test compareFields(dat, unpckd, skip=[:variableType])
@test compareFields(getStateKind(dat), getStateKind(unpckd))
@test isa(getStateKind(dat), ContinuousScalar)
@test isa(getStateKind(unpckd), ContinuousScalar)

##
end


@testset "test serialization of ManifoldKernelDensity" begin
##

# create a basic manifoldkerneldensity
mkd = manikde!(TranslationGroup(2), [randn(2) for _ in 1:100])

# convert up and down
st = pack(mkd) # TODO, PackedBelief
upk = unpack(st)

# and check the basics
@test isapprox( getPoints(mkd)[1], getPoints(upk)[1])
@test isapprox( getPoints(mkd)[end], getPoints(upk)[end])

@test mkd.manifold == upk.manifold
@test mkd._partial == upk._partial
@test mkd.infoPerCoord == upk.infoPerCoord

##
end



#
