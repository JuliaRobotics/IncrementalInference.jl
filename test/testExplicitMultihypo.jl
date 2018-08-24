
# Test multihypo computation assembly

using Distributions
using DocStringExtensions
using Base: Test

using IncrementalInference

@testset "test IncrementalInference.assembleHypothesesElements! without multihypothesis..." begin

s2_1_gt1 = 1:2
s2_1_gt2 = (1:20,Int[])
s2_1_gt3 = [(1,1:2); (2,Int[])]
s2_1_gt4 = Int[]
s2_1 = IncrementalInference.assembleHypothesesElements!(nothing, 20, 1, 2 )
@test sum(s2_1_gt1 .- s2_1[1]) == 0
@test sum(s2_1_gt2[1] .- s2_1[2][1]) == 0
@test sum(s2_1_gt2[2] .- s2_1[2][2]) == 0
@test s2_1_gt3[1][1] == s2_1[3][1][1]
@test sum(s2_1_gt3[1][2] .- s2_1[3][1][2]) == 0
@test s2_1_gt3[2][1] == s2_1[3][2][1]
@test sum(s2_1_gt3[2][2] .- s2_1[3][2][2]) == 0
@test sum(s2_1_gt4 .- s2_1[4])  == 0


s2_2_gt1 = 1:2
s2_2_gt2 = (1:20,Int[])
s2_2_gt3 = [(1,1:2); (2,Int[])]
s2_2_gt4 = Int[]

s2_2 = IncrementalInference.assembleHypothesesElements!(nothing, 20, 2, 2 )

@test sum(s2_2_gt1 .- s2_2[1]) == 0
@test sum(s2_2_gt2[1] .- s2_2[2][1]) == 0
@test sum(s2_2_gt2[2] .- s2_2[2][2]) == 0
@test s2_2_gt3[1][1] == s2_2[3][1][1]
@test sum(s2_2_gt3[1][2] .- s2_2[3][1][2]) == 0
@test s2_2_gt3[2][1] == s2_2[3][2][1]
@test sum(s2_2_gt3[2][2] .- s2_2[3][2][2]) == 0
@test sum(s2_2_gt4 .- s2_2[4])  == 0

end



@testset "test IncrementalInference.assembleHypothesesElements! with bi-modality..." begin

s3_1_gt1 = [1]
s3_1_gt2 = (0,5,5,20)
s3_1_gt3 = [(1,Int[]); (2,Int[1;2]); (3,Int[1;3])]
s3_1_gt4 = 20

s3_1 = IncrementalInference.assembleHypothesesElements!(Categorical([0.0;0.5;0.5]), 20, 1, 3 )

@test sum(s3_1_gt1 - s3_1[1]) == 0
@test sum(s3_1_gt2[1] .- s3_1[2][1]) == 0
@test length(s3_1[2][2]) > s3_1_gt2[2]
@test length(s3_1[2][3]) > s3_1_gt2[3]
@test length(s3_1[2][2]) + length(s3_1[2][3]) == s3_1_gt2[4]
@test s3_1_gt3[1][1] == s3_1[3][1][1]
@test s3_1_gt3[2][1] == s3_1[3][2][1]
@test s3_1_gt3[3][1] == s3_1[3][3][1]
@test sum(s3_1_gt3[1][2] .- s3_1[3][1][2]) == 0
@test sum(s3_1_gt3[2][2] .- s3_1[3][2][2]) == 0
@test sum(s3_1_gt3[3][2] .- s3_1[3][3][2]) == 0

@test sum(s3_1[4] .== 2) > s3_1_gt2[2]
@test sum(s3_1[4] .== 3) > s3_1_gt2[3]

@test sum( [1:20;][s3_1[4] .== 2] .== s3_1[2][2] ) == length(s3_1[2][2])
@test sum( [1:20;][s3_1[4] .== 3] .== s3_1[2][3] ) == length(s3_1[2][3])
@test length(s3_1[4]) == s3_1_gt4



s3_2_gt1 = [1]
s3_2_gt2 = (0,5,5,20)
s3_2_gt3 = [(1,Int[1;2]); (2,Int[2;3]); (3,Int[2;3])]
s3_2_gt4 = 20

s3_2 = IncrementalInference.assembleHypothesesElements!(Categorical([0.0;0.5;0.5]), 20, 2, 3 )

@test sum(s3_2_gt1 - s3_2[1]) == 0
@test sum(s3_2_gt2[1] .- s3_2[2][1]) == 0
@test length(s3_2[2][2]) > s3_2_gt2[2]
@test length(s3_2[2][3]) > s3_2_gt2[3]
@test length(s3_2[2][2]) + length(s3_2[2][3]) == s3_2_gt2[4]
@test s3_2_gt3[1][1] == s3_2[3][1][1]
@test s3_2_gt3[2][1] == s3_2[3][2][1]
@test s3_2_gt3[3][1] == s3_2[3][3][1]
@test sum(s3_2_gt3[1][2] .- s3_2[3][1][2]) == 0
@test sum(s3_2_gt3[2][2] .- s3_2[3][2][2]) == 0
@test sum(s3_2_gt3[3][2] .- s3_2[3][3][2]) == 0

@test sum(s3_2[4] .== 2) > s3_2_gt2[2]
@test sum(s3_2[4] .== 3) > s3_2_gt2[3]

@test sum( [1:20;][s3_2[4] .== 2] .== s3_2[2][2] ) == length(s3_2[2][2])
@test sum( [1:20;][s3_2[4] .== 3] .== s3_2[2][3] ) == length(s3_2[2][3])
@test length(s3_2[4]) == s3_2_gt4



s3_3_gt1 = [1]
s3_3_gt2 = (0,5,5,20)
s3_3_gt3 = [(1,Int[1;3]); (2,Int[2;3]); (3,Int[2;3])]
s3_3_gt4 = 20

s3_3 = IncrementalInference.assembleHypothesesElements!(Categorical([0.0;0.5;0.5]), 20, 3, 3 )

@test sum(s3_3_gt1 - s3_3[1]) == 0
@test sum(s3_3_gt2[1] .- s3_3[2][1]) == 0
@test length(s3_3[2][2]) > s3_3_gt2[2]
@test length(s3_3[2][3]) > s3_3_gt2[3]
@test length(s3_3[2][2]) + length(s3_3[2][3]) == s3_3_gt2[4]
@test s3_3_gt3[1][1] == s3_3[3][1][1]
@test s3_3_gt3[2][1] == s3_3[3][2][1]
@test s3_3_gt3[3][1] == s3_3[3][3][1]
@test sum(s3_3_gt3[1][2] .- s3_3[3][1][2]) == 0
@test sum(s3_3_gt3[2][2] .- s3_3[3][2][2]) == 0
@test sum(s3_3_gt3[3][2] .- s3_3[3][3][2]) == 0

@test sum(s3_3[4] .== 2) > s3_3_gt2[2]
@test sum(s3_3[4] .== 3) > s3_3_gt2[3]

@test sum( [1:20;][s3_3[4] .== 2] .== s3_3[2][2] ) == length(s3_3[2][2])
@test sum( [1:20;][s3_3[4] .== 3] .== s3_3[2][3] ) == length(s3_3[2][3])
@test length(s3_3[4]) == s3_3_gt4


end



@testset "test IncrementalInference.assembleHypothesesElements! with tri-modality..." begin


# s3_4_gt1 = [1]
# s3_4_gt2 = (0,5,5,20)
# s3_4_gt3 = [(1,Int[1;3]); (2,Int[2;3]); (3,Int[2;3])]
# s3_4_gt4 = 20

s4_2 = IncrementalInference.assembleHypothesesElements!(Categorical([0.0;0.3;0.3;0.4]), 20, 2, 4 )



end
