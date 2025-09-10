using IncrementalInference
using InteractiveUtils
using Test

##

@testset "test the basics" begin
##

fg = initfg()

addVariable!(fg, :x1, ContinuousScalar)
addVariable!(fg, :x2, ContinuousScalar)
addFactor!(fg, [:x1;:x2], LinearRelative(Normal()), graphinit=false)
addFactor!(fg, [:x2], Prior(Normal()), graphinit=false)

@test exists(fg, :x1)

@test !exists(fg, :l13)

##
end


@testset "test manikde! constructions on variableType" begin
##

pts = [randn(1) for _ in 1:100]
varT = LinearRelative(Normal(1.0))
manikde!(varT, pts)


DFG.@defStateType _TestManiKde SpecialEuclideanGroup(2; variant=:right) ArrayPartition([0;0.], [1 0; 0 1.])

# construct directly with ArrayPartition
pts = [ArrayPartition(randn(2), [1 0; 0 1.]) for _ in 1:100]
varT = _TestManiKde
manikde!(varT, pts)

# construct indirectly via tuple (expect users only, not meant for general use)
pts = [(randn(2), [1 0; 0 1.]) for _ in 1:100]
varT = _TestManiKde
manikde!(varT, pts)


##
end


@testset "test InteractiveUtilsExt" begin
##

IIF.listTypeTree(RelativeObservation)

IIF.getCurrentWorkspaceFactors()
IIF.getCurrentWorkspaceVariables()


##
end

#