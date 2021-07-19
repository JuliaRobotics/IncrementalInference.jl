
using IncrementalInference
using Test

##

@testset "test _evalFactorTemporary" begin
## test utility to build a temporary graph

fct = EuclidDistance(Normal(10,1))
T_pt_vec = [(ContinuousScalar,[0;]); (ContinuousScalar,[9.5;])]

##

dfg, _dfgfct = IIF._buildGraphByFactorAndTypes!(fct, T_pt_vec...)


## test  the evaluation of factor without

B = IIF._evalFactorTemporary!(EuclidDistance(Normal(10,1)), 2, ([10;],), T_pt_vec... );

@test_broken B isa Vector{Vector{Float64}}
@test isapprox( B[1], [10.0;], atol=1e-6)

##
end


@testset "test residual slack prerequisite for numerical factor gradients, Euclidean(1)" begin
##

fct = EuclidDistance(Normal(10,1))
measurement = ([10;],)
T_pt_args = [(ContinuousScalar,[0;]); (ContinuousScalar,[9.5;])]

##

slack_resid = calcFactorResidualTemporary(fct, measurement, T_pt_args...)

##

coord_1 = IIF._evalFactorTemporary!(fct, 1, measurement, T_pt_args..., _slack=slack_resid )
@test length(coord_1) == 1
@test isapprox( coord_1[1], [0.0], atol=1e-6)

coord_2 = IIF._evalFactorTemporary!(fct, 2, measurement, T_pt_args..., _slack=slack_resid )
@test length(coord_2) == 1
@test isapprox( coord_2[1], [9.5], atol=1e-6)

##

coord_1 = IIF._evalFactorTemporary!(fct, 1, measurement, T_pt_args... )
@test length(coord_1) == 1
@test isapprox( coord_1[1], [-0.5], atol=1e-6)

coord_2 = IIF._evalFactorTemporary!(fct, 2, measurement, T_pt_args... )
@test length(coord_2) == 1
@test isapprox( coord_2[1], [10.0], atol=1e-6)


##
end


@testset "test residual slack prerequisite for numerical factor gradients, Euclidean(2)" begin
##

pp = LinearRelative(MvNormal([10;0],[1 0; 0 1]))
measurement = ([10.0;0.0],)
T_pt_args = [(ContinuousEuclid{2},[0;0.0]); (ContinuousEuclid{2},[9.5;0])]

## test the building of factor graph to be correct

_fg,_ = IIF._buildGraphByFactorAndTypes!(pp, T_pt_args...);

@test length(getVal(_fg[:x1])) == 1
@test length(getVal(_fg[:x1])[1]) == 2

@test length(getVal(_fg[:x2])) == 1
@test length(getVal(_fg[:x2])[1]) == 2


##

_fg = initfg()

slack_resid = calcFactorResidualTemporary(fct, measurement, T_pt_args..., tfg=_fg)

@test length(getVal(_fg[:x1])) == 1
@test length(getVal(_fg[:x1])[1]) == 2

@test length(getVal(_fg[:x2])) == 1
@test length(getVal(_fg[:x2])[1]) == 2

##

tfg,_ = IIF._buildGraphByFactorAndTypes!(pp, T_pt_args...);

coord_1 = IIF._evalFactorTemporary!(fct, 1, measurement, T_pt_args..., _slack=slack_resid, tfg=tfg, newFactor=false )

##
@test length(coord_1) == 1
@test isapprox( coord_1[1], [0.0], atol=1e-6)

coord_2 = IIF._evalFactorTemporary!(fct, 2, measurement, T_pt_args..., _slack=slack_resid )
@test length(coord_2) == 1
@test isapprox( coord_2[1], [9.5], atol=1e-6)

##

coord_1 = IIF._evalFactorTemporary!(fct, 1, measurement, T_pt_args... )
@test length(coord_1) == 1
@test isapprox( coord_1[1], [-0.5], atol=1e-6)

coord_2 = IIF._evalFactorTemporary!(fct, 2, measurement, T_pt_args... )
@test length(coord_2) == 1
@test isapprox( coord_2[1], [10.0], atol=1e-6)


##
end


#