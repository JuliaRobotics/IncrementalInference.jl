using DistributedFactorGraphs
using IncrementalInference

using Test


##

mutable struct MutableLinearConditional{N, T <: SamplableBelief} <: AbstractRelativeFactor
    Z::T
    # timestamp::DateTime
end

function MutableLinearConditional{N}() where N
    newval = MvNormal(zeros(N), diagm(ones(N)))
    MutableLinearConditional{N,typeof(newval)}(newval)
end
MutableLinearConditional(n::Int=1) = MutableLinearConditional{n}()
MutableLinearConditional(nm::Distributions.ContinuousUnivariateDistribution) = MutableLinearConditional{1, typeof(nm)}(nm)
MutableLinearConditional(nm::MvNormal) = MutableLinearConditional{length(nm.μ), typeof(nm)}(nm)
MutableLinearConditional(nm::BallTreeDensity) = MutableLinearConditional{Ndim(nm), typeof(nm)}(nm)

getDimension(::Type{MutableLinearConditional{N,<:SamplableBelief}}) where {N} = N
getManifolds(::Type{MutableLinearConditional{N,<:SamplableBelief}}) where {N} = tuple([:Euclid for i in 1:N]...)

IIF.getSample(s::MutableLinearConditional, N::Int=1) = (reshape(rand(s.Z,N),:,N), )
function (s::MutableLinearConditional)(res::AbstractArray{<:Real},
                                userdata::FactorMetadata,
                                idx::Int,
                                meas::Tuple{<:AbstractArray{<:Real, 2}},
                                X1::AbstractArray{<:Real,2},
                                X2::AbstractArray{<:Real,2}  )
#
res[:] = meas[1][:,idx] - (X2[:,idx] - X1[:,idx])
nothing
end


IIF.getFactorMean(fct::MutableLinearConditional) = getFactorMean(fct.Z)
##


@testset "testing dead reckoning tether" begin

# start with an empty factor graph object
fg = initfg()

# Add the first pose :x0
addVariable!(fg, :x0, ContinuousScalar)

# Add at a fixed location PriorPose2 to pin :x0 to a starting location (10,10, pi/4)
addFactor!(fg, [:x0], Prior( Normal(0.0,0.1) ))

# Drive around in line
for i in 0:5
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, ContinuousScalar)
    pp = LinearConditional(Normal(1.0,0.1))
    addFactor!(fg, [psym;nsym], pp )
end

# Add landmarks with Bearing range measurements
addVariable!(fg, :l1, ContinuousScalar, tags=[:LANDMARK])
p2br = LinearConditional(Normal(1.0,0.1))
addFactor!(fg, [:x0; :l1], p2br )


## async solving with dead-reckon branch

addVariable!(fg, :deadreckon_x0, ContinuousScalar, solvable=0)

drec = MutableLinearConditional(Normal(0.0,0.1))

addFactor!(fg, [:x0; :deadreckon_x0], drec, solvable=0)

#
@test length(map( x->x.label, getVariables(fg, solvable=1))) == 8
@test length(map( x->x.label, getVariables(fg, solvable=0))) == 9
#
# # make sure
@test length(getEliminationOrder(fg, solvable=1)) == 8
# check default
@test length(getEliminationOrder(fg)) == 8

# default check
vo = getEliminationOrder(fg)
@test length(vo) == 8


tree = buildTreeFromOrdering!(fg,vo)


tree2, smt, hists = solveTree!(fg);

@test !isInitialized(fg, :deadreckon_x0)


val = accumulateFactorMeans(fg, [:x0deadreckon_x0f1])

@test norm(val - calcVariablePPE(fg, :x0).suggested) < 1e-8


end
