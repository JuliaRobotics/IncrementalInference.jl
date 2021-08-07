

using Test
using IncrementalInference

##

import IncrementalInference: getSample

mutable struct PartialDim2{T} <: AbstractPrior
  Z::T
  partial::Tuple{Int}
end

PartialDim2(Z::D) where {D <: IIF.SamplableBelief} = PartialDim2(Z, (2,))

getSample(cfo::CalcFactor{<:PartialDim2}, N::Int=1) = ([rand(cfo.factor.Z, 1)[:] for _ in 1:N],)


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

predictbelief(fg, :x0, [:x0f1;])
@test true

##

end