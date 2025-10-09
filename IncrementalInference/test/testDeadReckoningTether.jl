using DistributedFactorGraphs
using IncrementalInference

using Test


##

mutable struct MutableLinearRelative{N, T <: SamplableBelief} <: RelativeObservation
    Z::T
    # timestamp::DateTime
end

function MutableLinearRelative{N}() where N
    newval = MvNormal(zeros(N), diagm(ones(N)))
    MutableLinearRelative{N,typeof(newval)}(newval)
end
MutableLinearRelative(n::Int=1) = MutableLinearRelative{n}()
MutableLinearRelative(nm::Distributions.ContinuousUnivariateDistribution) = MutableLinearRelative{1, typeof(nm)}(nm)
MutableLinearRelative(nm::MvNormal) = MutableLinearRelative{length(nm.μ), typeof(nm)}(nm)
MutableLinearRelative(nm::ManifoldKernelDensity) = MutableLinearRelative{Ndim(nm), typeof(nm)}(nm)

DFG.getDimension(::Type{MutableLinearRelative{N,<:SamplableBelief}}) where {N} = N
DFG.getManifold(::MutableLinearRelative{N}) where N = TranslationGroup(N)


function IIF.getSample(cf::CalcFactor{<:MutableLinearRelative})
    return rand(cf.factor.Z, 1)
end

function (s::CalcFactor{<:MutableLinearRelative})(  meas,
                                                    X1,
                                                    X2  )
#
    return meas .- (X2 .- X1)
end

##

@testset "testing dead reckoning tether" begin

# test error message and then define method for MutableLinearRelative
# @test_throws ErrorException getFactorMean(MutableLinearRelative(Normal(0.0,0.1)))

# IIF.getFactorMean(fct::MutableLinearRelative) = getFactorMean(fct.Z)
##

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
    pp = LinearRelative(Normal(1.0,0.1))
    addFactor!(fg, [psym;nsym], pp )
end

# Add landmarks with Bearing range measurements
addVariable!(fg, :l1, ContinuousScalar, tags=[:LANDMARK])
p2br = LinearRelative(Normal(1.0,0.1))
addFactor!(fg, [:x0; :l1], p2br )


## async solving with dead-reckon branch

addVariable!(fg, :deadreckon_x0, ContinuousScalar, solvable=0)

drec = MutableLinearRelative(Normal(0.0,0.1))

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

tree2 = solveTree!(fg);

@test !isInitialized(fg, :deadreckon_x0)

val = accumulateFactorMeans(fg, [:x0deadreckon_x0f1])

# must fix return type stability
fval = float(val...)
@test isapprox(fval, calcMeanMaxSuggested(fg, :x0).suggested[1], atol=1e-4 )

#TODO improve test
rebaseFactorVariable!(fg, :x0deadreckon_x0f1, [:x1; :deadreckon_x0])
@test issetequal(ls2(fg, :x0), [:x1, :l1])
@test issetequal(ls2(fg, :x1), [:x0, :deadreckon_x0, :x2])

##

end
