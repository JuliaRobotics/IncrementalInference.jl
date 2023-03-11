
using IncrementalInference
using Test

##

@testset "test _evalFactorTemporary" begin
## test utility to build a temporary graph

fct = EuclidDistance(Normal(10,1))
varTypes = (ContinuousScalar,ContinuousScalar); 
varPts = ([0;],[9.5;])

##

dfg, _dfgfct = IIF._buildGraphByFactorAndTypes!(fct, varTypes, varPts)

@test length(intersect(ls(dfg), [:x1; :x2])) == 2
@test lsf(dfg) == [:x1x2f1;]

## test  the evaluation of factor without

B = IIF._evalFactorTemporary!(EuclidDistance(Normal(10,1)), varTypes, 2, [[10;]], varPts );

@test B isa Vector{Vector{Float64}}
@test isapprox( B[1], [10.0;], atol=1e-6)

##
end


@testset "test residual slack prerequisite for numerical factor gradients, Euclidean(1)" begin
##

fct = EuclidDistance(Normal(10,1))
measurement = [[10;]]
varTypes = (ContinuousScalar,ContinuousScalar)
pts = ([0;],[9.5;])

##

slack_resid = calcFactorResidualTemporary(fct, varTypes, measurement[1], pts)

##

coord_1 = IIF._evalFactorTemporary!(fct, varTypes, 1, measurement, pts, _slack=slack_resid )
@test length(coord_1) == 1
@test isapprox( coord_1[1], [0.0], atol=1e-6)

coord_2 = IIF._evalFactorTemporary!(fct, varTypes, 2, measurement, pts, _slack=slack_resid )
@test length(coord_2) == 1
@test isapprox( coord_2[1], [9.5], atol=1e-6)

##

coord_1 = IIF._evalFactorTemporary!(fct, varTypes, 1, measurement, pts )
@test length(coord_1) == 1
@test isapprox( coord_1[1], [-0.5], atol=1e-6)

coord_2 = IIF._evalFactorTemporary!(fct, varTypes, 2, measurement, pts )
@test length(coord_2) == 1
@test isapprox( coord_2[1], [10.0], atol=1e-6)


##
end


@testset "test residual slack prerequisite for numerical factor gradients, Euclidean(2)" begin
##

fct = LinearRelative(MvNormal([10;0],[1 0; 0 1]))
measurement = [[10.0;0.0]]
varTypes = (ContinuousEuclid{2},ContinuousEuclid{2})
pts = ([0;0.0] ,[9.5;0])

## test the building of factor graph to be correct

_fg,_ = IIF._buildGraphByFactorAndTypes!(fct, varTypes, pts);

@test length(getVal(_fg[:x1])) == 1
@test length(getVal(_fg[:x1])[1]) == 2

@test length(getVal(_fg[:x2])) == 1
@test length(getVal(_fg[:x2])[1]) == 2


##

_fg = initfg()

slack_resid = calcFactorResidualTemporary(fct, varTypes, measurement[1], pts, tfg=_fg)

@test length(getVal(_fg[:x1])) == 1
@test length(getVal(_fg[:x1])[1]) == 2

@test length(getVal(_fg[:x2])) == 1
@test length(getVal(_fg[:x2])[1]) == 2

## Manually provide a common temp graph and force no factor and same variables via keywords

tfg,_ = IIF._buildGraphByFactorAndTypes!(fct, varTypes, pts);
coord_1 = IIF._evalFactorTemporary!(fct, varTypes, 1, measurement, pts, _slack=slack_resid, tfg=tfg, newFactor=false, currNumber=0 )

##

@test length(coord_1) == 1
@test length(coord_1[1]) == 2

@test isapprox( coord_1[1], [0;0.0], atol=1e-6)

coord_2 = IIF._evalFactorTemporary!(fct, varTypes, 2, measurement, pts, _slack=slack_resid )
@test length(coord_2) == 1
@test length(coord_2[1]) == 2

@test isapprox( coord_2[1], [9.5; 0], atol=1e-6)

## Repeat the same test but allow _evalFactorTemporary to self construct internal temporary graph

coord_1 = IIF._evalFactorTemporary!(fct, varTypes, 1, measurement, pts )

##
@test length(coord_1) == 1
@test length(coord_1[1]) == 2

@test isapprox( coord_1[1], [-0.5;0], atol=1e-6)

coord_2 = IIF._evalFactorTemporary!(fct, varTypes, 2, measurement, pts )
@test length(coord_2) == 1
@test length(coord_2[1]) == 2

@test isapprox( coord_2[1], [10.0;0], atol=1e-6)

##
end


@testset "Enable SolverParams.attemptGradients" begin
##

fg = generateGraph_LineStep(4;
  solverParams = SolverParams(;
    attemptGradients=true
  )
)

##
end

#