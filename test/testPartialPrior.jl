

using Test
using IncrementalInference
using Manifolds
using DistributedFactorGraphs
import IncrementalInference: getSample, getManifold

##

mutable struct PartialDim2{T} <: AbstractPrior
  Z::T
  partial::Tuple{Int}
end

PartialDim2(Z::D) where {D <: IIF.SamplableBelief} = PartialDim2(Z, (2,))

getManifold(pd::PartialDim2) = getManifoldPartial(TranslationGroup(2), [pd.partial...])[1]

getSample(cfo::CalcFactor{<:PartialDim2}) = [rand(cfo.factor.Z);]


##

@testset "Test partial dimensions on prior are correct" begin
##

fg = initfg()

v0 = addVariable!(fg, :x0, ContinuousEuclid{2})

f0 = addFactor!(fg, [:x0], PartialDim2(Normal()))

@test IIF._getDimensionsPartial(f0) == [2]

bel, infd = propagateBelief(fg, v0, [f0;])

@test isPartial(bel)

##

dens = Vector{ManifoldKernelDensity}()
IIF.proposalbeliefs!(fg, :x0, [f0], dens)

pts = getPoints(dens[1], false)

##

propagateBelief(fg, :x0, [:x0f1;])
@test true

##
end



@testset "test propagateBelief returning a partial" begin
##

fg = initfg()

v0 = addVariable!(fg, :x0, ContinuousEuclid{2})

pts = [randn(1) for _ in 1:1000];
mkd = manikde!(TranslationGroup(1), pts, bw=[0.1;])
pp = PartialPrior(ContinuousEuclid{2}, mkd, (2,))
f0 = addFactor!(fg, [:x0;], pp, graphinit=false)

##

bel, infd = propagateBelief(fg, v0, [f0;])
@test isPartial(bel)

##

doautoinit!(fg, :x0)

# check the number of points in the graph value store
@show getSolverParams(fg).N
@test length(getPoints(getBelief(fg, :x0))) == getSolverParams(fg).N
@info "PassThrough factors currently work different and will pass the full N=1000 count through to the graph."

##
end
