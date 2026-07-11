
# Test multihypo computation assembly

using Test
using IncrementalInference


@testset "test IncrementalInference._prepareHypoRecipe! with only nullhypothesis..." begin

##

# n2_1 == (certainidx, allelements, activehypo, mhidx)
n2_1_gt1 = 1:2
n2_1_gt2_ = (3,3,0)
n2_1_gt3 = [(0,Int[1;]); (1,1:2); (2,Int[])]
n2_1_gt4_ = 20
n2_1 = IncrementalInference._prepareHypoRecipe!(nothing, 20, 1, 2, ones(Bool, 2), 0.5 )

@test sum([n2_1_gt1;] - [n2_1.certainidx;]) == 0
@test length(n2_1.allelements[1]) > n2_1_gt2_[1]
@test length(n2_1.allelements[2]) > n2_1_gt2_[2]
@test length(n2_1.allelements[3]) == n2_1_gt2_[3]
@test length(n2_1.allelements[1]) + length(n2_1.allelements[2]) == n2_1_gt4_
@test n2_1_gt3[1][1] == n2_1.activehypo[1][1]
@test n2_1_gt3[2][1] == n2_1.activehypo[2][1]
@test n2_1_gt3[3][1] == n2_1.activehypo[3][1]
@test sum(n2_1_gt3[1][2] .- n2_1.activehypo[1][2]) == 0
@test sum([n2_1_gt3[2][2];] .- [n2_1.activehypo[2][2];]) == 0
@test sum(n2_1_gt3[3][2] .- n2_1.activehypo[3][2]) == 0
@test sum(n2_1.mhidx .== 0) > n2_1_gt2_[1]
@test sum(n2_1.mhidx .== 1) > n2_1_gt2_[1]
@test sum( [1:n2_1_gt4_;][n2_1.mhidx .== 0] .== n2_1.allelements[1] ) == length(n2_1.allelements[1])
@test sum( [1:n2_1_gt4_;][n2_1.mhidx .== 1] .== n2_1.allelements[2] ) == length(n2_1.allelements[2])
@test length(n2_1.mhidx) == n2_1_gt4_



# n2_1 == (certainidx, allelements, activehypo, mhidx)
n2_1_gt1 = 1:2
n2_1_gt2_ = (3,3,0)
n2_1_gt3 = [(0,Int[2;]); (1,1:2); (2,Int[])]
n2_1_gt4_ = 20
n2_1 = IncrementalInference._prepareHypoRecipe!(nothing, 20, 2, 2, ones(Bool, 2), 0.5 )

@test sum([n2_1_gt1;] - [n2_1.certainidx;]) == 0
@test length(n2_1.allelements[1]) > n2_1_gt2_[1]
@test length(n2_1.allelements[2]) > n2_1_gt2_[2]
@test length(n2_1.allelements[3]) == n2_1_gt2_[3]
@test length(n2_1.allelements[1]) + length(n2_1.allelements[2]) == n2_1_gt4_
@test n2_1_gt3[1][1] == n2_1.activehypo[1][1]
@test n2_1_gt3[2][1] == n2_1.activehypo[2][1]
@test n2_1_gt3[3][1] == n2_1.activehypo[3][1]
@test sum(n2_1_gt3[1][2] .- n2_1.activehypo[1][2]) == 0
@test sum([n2_1_gt3[2][2];] .- [n2_1.activehypo[2][2];]) == 0
@test sum(n2_1_gt3[3][2] .- n2_1.activehypo[3][2]) == 0
@test sum(n2_1.mhidx .== 0) > n2_1_gt2_[1]
@test sum(n2_1.mhidx .== 1) > n2_1_gt2_[1]
@test sum( [1:n2_1_gt4_;][n2_1.mhidx .== 0] .== n2_1.allelements[1] ) == length(n2_1.allelements[1])
@test sum( [1:n2_1_gt4_;][n2_1.mhidx .== 1] .== n2_1.allelements[2] ) == length(n2_1.allelements[2])
@test length(n2_1.mhidx) == n2_1_gt4_

##

end


@testset "test IncrementalInference._prepareHypoRecipe! without multihypothesis..." begin

##

# certainidx = 1 ## ??
# sfidx=1, mhidx=0:  ah = [1;]
# sfidx=1, mhidx=1:  ah = [1;2]

# s2_1 == (certainidx, allelements, activehypo, mhidx)
s2_1_gt1 = 1:2
s2_1_gt2 = (Int[],1:20,Int[])
s2_1_gt3 = [(0,Int[1;]); (1,1:2); (2,Int[])]
s2_1_gt4 = ones(20)
s2_1 = IncrementalInference._prepareHypoRecipe!(nothing, 20, 1, 2 )
@test sum([s2_1_gt1;] .- [s2_1.certainidx;]) == 0
@test sum( s2_1_gt2[1] .- s2_1.allelements[1]) == 0
@test sum([s2_1_gt2[2];] .- [s2_1.allelements[2];]) == 0
@test sum( s2_1_gt2[3] .- s2_1.allelements[3]) == 0
@test s2_1_gt3[1][1] == s2_1.activehypo[1][1]
@test sum(s2_1_gt3[1][2] .- s2_1.activehypo[1][2]) == 0
@test s2_1_gt3[2][1] == s2_1.activehypo[2][1]
@test sum([s2_1_gt3[2][2];] .- [s2_1.activehypo[2][2];]) == 0
@test s2_1_gt3[3][1] == s2_1.activehypo[3][1]
@test sum(s2_1_gt3[3][2] .- s2_1.activehypo[3][2]) == 0
@test sum(s2_1_gt4 .- s2_1.mhidx)  == 0


s2_2_gt1 = 1:2
s2_2_gt2 = (Int[],1:20,Int[])
s2_2_gt3 = [(0,Int[2;]); (1,1:2); (2,Int[])]
s2_2_gt4 = ones(20) # Int[]

s2_2 = IncrementalInference._prepareHypoRecipe!(nothing, 20, 2, 2 )

@test sum([s2_2_gt1;] .- [s2_2.certainidx;]) == 0
@test sum([s2_2_gt2[1];] .- [s2_2.allelements[1];]) == 0
@test sum(s2_2_gt2[2] .- s2_2.allelements[2]) == 0
@test sum(s2_2_gt2[3] .- s2_2.allelements[3]) == 0
@test      s2_2_gt3[1][1]   ==  s2_2.activehypo[1][1]
@test sum([s2_2_gt3[1][2];] .- [s2_2.activehypo[1][2];]) == 0
@test s2_2_gt3[2][1] == s2_2.activehypo[2][1]
@test sum([s2_2_gt3[2][2];] .- [s2_2.activehypo[2][2];]) == 0
@test s2_2_gt3[3][1] == s2_2.activehypo[3][1]
@test sum([s2_2_gt3[3][2];] .- [s2_2.activehypo[3][2];]) == 0
@test sum(s2_2_gt4 .- s2_2.mhidx)  == 0

##

end



@testset "_prepareHypoRecipe! with bi-modality (certain variable)" begin

##

# certainidx = 1
# sfidx=1, mhidx=1:  ah = []
# sfidx=1, mhidx=2:  ah = [1;2]
# sfidx=1, mhidx=3:  ah = [1;3]
s3_1_gt1 = [1]
s3_1_gt2 = (0,3,3,40)
s3_1_gt3 = [(1,Int[]); (2,Int[1;2]); (3,Int[1;3])]
s3_1_gt4 = 40

s3_1 = IncrementalInference._prepareHypoRecipe!(Categorical([0.0;0.5;0.5]), 40, 1, 3)

@test sum([s3_1_gt1;] - [s3_1.certainidx;]) == 0
@test sum([s3_1_gt2[1];] .- [s3_1.allelements[1];]) == 0
@test length(s3_1.allelements[2]) > s3_1_gt2[2]
@test length(s3_1.allelements[3]) > s3_1_gt2[3]
@test length(s3_1.allelements[2]) + length(s3_1.allelements[3]) == s3_1_gt2[4]
@test s3_1_gt3[1][1] == s3_1.activehypo[1][1]
@test s3_1_gt3[2][1] == s3_1.activehypo[2][1]
@test s3_1_gt3[3][1] == s3_1.activehypo[3][1]
@test sum(s3_1_gt3[1][2] .- s3_1.activehypo[1][2]) == 0
@test sum(s3_1_gt3[2][2] .- s3_1.activehypo[2][2]) == 0
@test sum(s3_1_gt3[3][2] .- s3_1.activehypo[3][2]) == 0

@test sum(s3_1.mhidx .== 2) > s3_1_gt2[2]
@test sum(s3_1.mhidx .== 3) > s3_1_gt2[3]

@test sum( [1:40;][s3_1.mhidx .== 2] .== s3_1.allelements[2] ) == length(s3_1.allelements[2])
@test sum( [1:40;][s3_1.mhidx .== 3] .== s3_1.allelements[3] ) == length(s3_1.allelements[3])
@test length(s3_1.mhidx) == s3_1_gt4

##

end


@testset "_prepareHypoRecipe! with bi-modality (fractional variable 1/2)" begin

##

# certainidx = 1
# sfidx=2, mhidx=1:  ah = [1;2]
# sfidx=2, mhidx=2:  ah = [2;3]
# sfidx=2, mhidx=3:  [2;3], 2 should take a value from 3
s3_2_gt1 = [1]
s3_2_gt2 = (0,3,3,40)
s3_2_gt3 = [(0, Int[2]); (1,Int[1;2]); (2,Int[1;2]); (3,Int[2;3])]
s3_2_gt4 = 40

s3_2 = IncrementalInference._prepareHypoRecipe!(Categorical([0.0;0.5;0.5]), 40, 2, 3 )

@test sum(s3_2_gt1 - s3_2.certainidx) == 0
@test sum(s3_2_gt2[1] .- s3_2.allelements[2]) == 0
@test length(s3_2.allelements[1]) > 0.5*s3_2_gt2[2] # reuse test reference for bad-init nullhypo case
@test length(s3_2.allelements[2]) == 0
@test length(s3_2.allelements[3]) > s3_2_gt2[2]
@test length(s3_2.allelements[4]) > s3_2_gt2[3]
@test length(s3_2.allelements[1]) + length(s3_2.allelements[3]) + length(s3_2.allelements[4]) == s3_2_gt2[4]
@test s3_2_gt3[1][1] == s3_2.activehypo[1][1]
@test s3_2_gt3[2][1] == s3_2.activehypo[2][1]
@test s3_2_gt3[3][1] == s3_2.activehypo[3][1]
@test sum(s3_2_gt3[1][2] .- s3_2.activehypo[1][2]) == 0
@test sum(s3_2_gt3[2][2] .- s3_2.activehypo[2][2]) == 0
@test sum(s3_2_gt3[3][2] .- s3_2.activehypo[3][2]) == 0

@test sum(s3_2.mhidx .== 2) > s3_2_gt2[2]
@test sum(s3_2.mhidx .== 3) > s3_2_gt2[3]

@test sum( [1:40;][s3_2.mhidx .== 0] .== s3_2.allelements[1] ) == length(s3_2.allelements[1])
@test 0 == length(s3_2.allelements[2])
@test sum( [1:40;][s3_2.mhidx .== 2] .== s3_2.allelements[3] ) == length(s3_2.allelements[3])
@test sum( [1:40;][s3_2.mhidx .== 3] .== s3_2.allelements[4] ) == length(s3_2.allelements[4])
@test length(s3_2.mhidx) == s3_2_gt4

##

end

@testset "_prepareHypoRecipe! with bi-modality (fractional variable 2/2)" begin

##

# certainidx = 1
# sfidx=3, mhidx=1:  ah = [1;3]
# sfidx=3, mhidx=2:  [2:3], 3 should take a value from 2
# sfidx=3, mhidx=3:  ah = [1;3]
s3_3_gt1 = [1]
s3_3_gt2 = (0,3,3,40)
s3_3_gt3 = [(0, Int[3]); (1,Int[1;3]); (2,Int[2;3]); (3,Int[1;3])]
s3_3_gt4 = 40

s3_3 = IncrementalInference._prepareHypoRecipe!(Categorical([0.0;0.5;0.5]), 40, 3, 3 )

@test sum(s3_3_gt1 - s3_3.certainidx) == 0
@test sum(s3_3_gt2[1] .- s3_3.allelements[2]) == 0
@test length(s3_3.allelements[1]) > 0.5*s3_3_gt2[2]
@test length(s3_3.allelements[2]) == 0
@test length(s3_3.allelements[3]) > s3_3_gt2[2]
@test length(s3_3.allelements[4]) > s3_3_gt2[3]
@test length(s3_3.allelements[1]) + length(s3_3.allelements[3]) + length(s3_3.allelements[4]) == s3_3_gt2[4]
@test s3_3_gt3[1][1] == s3_3.activehypo[1][1]
@test s3_3_gt3[2][1] == s3_3.activehypo[2][1]
@test s3_3_gt3[3][1] == s3_3.activehypo[3][1]
@test sum(s3_3_gt3[1][2] .- s3_3.activehypo[1][2]) == 0
@test sum(s3_3_gt3[2][2] .- s3_3.activehypo[2][2]) == 0
@test sum(s3_3_gt3[3][2] .- s3_3.activehypo[3][2]) == 0

@test sum(s3_3.mhidx .== 2) > s3_3_gt2[2]
@test sum(s3_3.mhidx .== 3) > s3_3_gt2[3]

@test sum( [1:40;][s3_3.mhidx .== 0] .== s3_3.allelements[1] ) == length(s3_3.allelements[1])
@test 0 == length(s3_3.allelements[2])
@test sum( [1:40;][s3_3.mhidx .== 2] .== s3_3.allelements[3] ) == length(s3_3.allelements[3])
@test sum( [1:40;][s3_3.mhidx .== 3] .== s3_3.allelements[4] ) == length(s3_3.allelements[4])
@test length(s3_3.mhidx) == s3_3_gt4

##

end



# @testset "test IncrementalInference._prepareHypoRecipe! with bi-modality backwards permutation..." begin

# certainidx = 1
# sfidx=1, mhidx=1:  ah = []
# sfidx=1, mhidx=2:  ah = [1;2]
# sfidx=1, mhidx=3:  ah = [1;3]
# s3_1_gt1 = [1]
# s3_1_gt2 = (0,3,3,20)
# s3_1_gt3 = [(1,Int[]); (2,Int[1;2]); (3,Int[1;3])]
# s3_1_gt4 = 20
#
# s3_1 = IncrementalInference._prepareHypoRecipe!(Categorical([0.5;0.5;0.0]), 20, 1, 3)

# @test sum(s3_1_gt1 - s3_1.certainidx) == 0
# @test sum(s3_1_gt2[1] .- s3_1.allelements[1]) == 0
# @test length(s3_1.allelements[2]) > s3_1_gt2[2]
# @test length(s3_1.allelements[3]) > s3_1_gt2[3]
# @test length(s3_1.allelements[2]) + length(s3_1.allelements[3]) == s3_1_gt2[4]
# @test s3_1_gt3[1][1] == s3_1.activehypo[1][1]
# @test s3_1_gt3[2][1] == s3_1.activehypo[2][1]
# @test s3_1_gt3[3][1] == s3_1.activehypo[3][1]
# @test sum(s3_1_gt3[1][2] .- s3_1.activehypo[1][2]) == 0
# @test sum(s3_1_gt3[2][2] .- s3_1.activehypo[2][2]) == 0
# @test sum(s3_1_gt3[3][2] .- s3_1.activehypo[3][2]) == 0
#
# @test sum(s3_1.mhidx .== 2) > s3_1_gt2[2]
# @test sum(s3_1.mhidx .== 3) > s3_1_gt2[3]
#
# @test sum( [1:20;][s3_1.mhidx .== 2] .== s3_1.allelements[2] ) == length(s3_1.allelements[2])
# @test sum( [1:20;][s3_1.mhidx .== 3] .== s3_1.allelements[3] ) == length(s3_1.allelements[3])
# @test length(s3_1.mhidx) == s3_1_gt4
#


# end


@testset "_prepareHypoRecipe! with tri-modality (certain variable)" begin

##

N = 50
s4_1_gt1 = [1]
s4_1_gt2 = (0,3,3,3,N)
s4_1_gt3 = [(1,Int[]); (2,Int[1;2]); (3,Int[1;3]); (4,Int[1;4])]
s4_1_gt4 = N

s4_1 = IncrementalInference._prepareHypoRecipe!(Categorical([0.0;0.33;0.33;0.34]), N, 1, 4 )

@test sum(s4_1_gt1 - s4_1.certainidx) == 0
@test sum(s4_1_gt2[1] .- s4_1.allelements[1]) == 0
@test length(s4_1.allelements[2]) > s4_1_gt2[2]
@test length(s4_1.allelements[3]) > s4_1_gt2[3]
@test length(s4_1.allelements[4]) > s4_1_gt2[4]
@test length(s4_1.allelements[2]) + length(s4_1.allelements[3]) + length(s4_1.allelements[4]) == s4_1_gt2[5]

@test s4_1_gt3[1][1] == s4_1.activehypo[1][1]
@test s4_1_gt3[2][1] == s4_1.activehypo[2][1]
@test s4_1_gt3[3][1] == s4_1.activehypo[3][1]
@test s4_1_gt3[4][1] == s4_1.activehypo[4][1]
@test sum(s4_1_gt3[1][2] .- s4_1.activehypo[1][2]) == 0
@test sum(s4_1_gt3[2][2] .- s4_1.activehypo[2][2]) == 0
@test sum(s4_1_gt3[3][2] .- s4_1.activehypo[3][2]) == 0
@test sum(s4_1_gt3[4][2] .- s4_1.activehypo[4][2]) == 0

@test sum(s4_1.mhidx .== 2) > s4_1_gt2[2]
@test sum(s4_1.mhidx .== 3) > s4_1_gt2[3]
@test sum(s4_1.mhidx .== 4) > s4_1_gt2[4]

@test sum( [1:N;][s4_1.mhidx .== 2] .== s4_1.allelements[2] ) == length(s4_1.allelements[2])
@test sum( [1:N;][s4_1.mhidx .== 3] .== s4_1.allelements[3] ) == length(s4_1.allelements[3])
@test sum( [1:N;][s4_1.mhidx .== 4] .== s4_1.allelements[4] ) == length(s4_1.allelements[4])
@test length(s4_1.mhidx) == s4_1_gt4

##


end


@testset "_prepareHypoRecipe! with tri-modality (fractional variable 1/3)" begin

## solve for fractional variable in trinary case

N = 70
s4_2_gt1 = [1]
s4_2_gt2 = (0,3,3,3,N)
s4_2_gt3 = [(0,Int[2]); (1,Int[1;2]); (2,Int[1;2]); (3,Int[2;3;4]); (4,Int[2;3;4])]
s4_2_gt4 = N

s4_2 = IncrementalInference._prepareHypoRecipe!(Categorical([0.0;0.33;0.33;0.34]), N, 2, 4 )

@test sum(s4_2_gt1 - s4_2.certainidx) == 0
@test length(s4_2.allelements[1]) > 0.5*s4_2_gt2[2]
@test sum(s4_2_gt2[2] .- s4_2.allelements[2]) == 0
@test length(s4_2.allelements[3]) > s4_2_gt2[2]
@test length(s4_2.allelements[4]) > s4_2_gt2[3]
@test length(s4_2.allelements[5]) > s4_2_gt2[4]
@test length(s4_2.allelements[1]) + length(s4_2.allelements[3]) + length(s4_2.allelements[4]) + length(s4_2.allelements[5]) == s4_2_gt2[5]

@test s4_2_gt3[1][1] == s4_2.activehypo[1][1]
@test s4_2_gt3[2][1] == s4_2.activehypo[2][1]
@test s4_2_gt3[3][1] == s4_2.activehypo[3][1]
@test s4_2_gt3[4][1] == s4_2.activehypo[4][1]
@test sum(s4_2_gt3[1][2] .- s4_2.activehypo[1][2]) == 0
@test sum(s4_2_gt3[2][2] .- s4_2.activehypo[2][2]) == 0
@test sum(s4_2_gt3[3][2] .- s4_2.activehypo[3][2]) == 0
@test sum(s4_2_gt3[4][2] .- s4_2.activehypo[4][2]) == 0

@test sum(s4_2.mhidx .== 0) > s4_2_gt2[2]
@test sum(s4_2.mhidx .== 2) > s4_2_gt2[2]
@test sum(s4_2.mhidx .== 3) > s4_2_gt2[3]
@test sum(s4_2.mhidx .== 4) > s4_2_gt2[4]

@test sum( [1:N;][s4_2.mhidx .== 0] .== s4_2.allelements[1] ) == length(s4_2.allelements[1])
@test 0 == length(s4_2.allelements[2])
@test sum( [1:N;][s4_2.mhidx .== 2] .== s4_2.allelements[3] ) == length(s4_2.allelements[3])
@test sum( [1:N;][s4_2.mhidx .== 3] .== s4_2.allelements[4] ) == length(s4_2.allelements[4])
@test sum( [1:N;][s4_2.mhidx .== 4] .== s4_2.allelements[5] ) == length(s4_2.allelements[5])
@test length(s4_2.mhidx) == s4_2_gt4

##

end



@testset "_prepareHypoRecipe! with tri-modality (fractional variable 2/3)" begin

##

N = 70
s4_3_gt1 = [1]
s4_3_gt2 = (0,3,3,3,N)
s4_3_gt3 = [(0,Int[3]); (1,Int[1;3]); (2,Int[2;3;4]); (3,Int[1;3]); (4,Int[2;3;4])]
s4_3_gt4 = N

s4_3 = IncrementalInference._prepareHypoRecipe!(Categorical([0.0;0.33;0.33;0.34]), N, 3, 4 )

@test sum(s4_3_gt1 - s4_3.certainidx) == 0
@test length(s4_3.allelements[1]) > 0.5*s4_3_gt2[2]
@test sum(s4_3_gt2[2] .- s4_3.allelements[2]) == 0
@test length(s4_3.allelements[3]) > s4_3_gt2[2]
@test length(s4_3.allelements[4]) > s4_3_gt2[3]
@test length(s4_3.allelements[5]) > s4_3_gt2[4]
@test length(s4_3.allelements[1]) + length(s4_3.allelements[3]) + length(s4_3.allelements[4]) + length(s4_3.allelements[5]) == s4_3_gt2[5]

@test s4_3_gt3[1][1] == s4_3.activehypo[1][1]
@test s4_3_gt3[2][1] == s4_3.activehypo[2][1]
@test s4_3_gt3[3][1] == s4_3.activehypo[3][1]
@test s4_3_gt3[4][1] == s4_3.activehypo[4][1]
@test sum(s4_3_gt3[1][2] .- s4_3.activehypo[1][2]) == 0
@test sum(s4_3_gt3[2][2] .- s4_3.activehypo[2][2]) == 0
@test sum(s4_3_gt3[3][2] .- s4_3.activehypo[3][2]) == 0
@test sum(s4_3_gt3[4][2] .- s4_3.activehypo[4][2]) == 0

@test sum(s4_3.mhidx .== 0) > s4_3_gt2[2]
@test sum(s4_3.mhidx .== 2) > s4_3_gt2[2]
@test sum(s4_3.mhidx .== 3) > s4_3_gt2[3]
@test sum(s4_3.mhidx .== 4) > s4_3_gt2[4]

@test sum( [1:N;][s4_3.mhidx .== 0] .== s4_3.allelements[1] ) == length(s4_3.allelements[1])
@test 0 == length(s4_3.allelements[2])
@test sum( [1:N;][s4_3.mhidx .== 2] .== s4_3.allelements[3] ) == length(s4_3.allelements[3])
@test sum( [1:N;][s4_3.mhidx .== 3] .== s4_3.allelements[4] ) == length(s4_3.allelements[4])
@test sum( [1:N;][s4_3.mhidx .== 4] .== s4_3.allelements[5] ) == length(s4_3.allelements[5])
@test length(s4_3.mhidx) == s4_3_gt4

##

end

@testset "_prepareHypoRecipe! with tri-modality (fractional variable 3/3)" begin

##

N = 70
s4_4_gt1 = [1]
s4_4_gt2 = (0,3,3,3,N)
s4_4_gt3 = [(0,Int[4]); (1,Int[1;4]); (2,Int[2;3;4]); (3,Int[2;3;4]); (4,Int[1;4])]
s4_4_gt4 = N

s4_4 = IncrementalInference._prepareHypoRecipe!(Categorical([0.0;0.33;0.33;0.34]), N, 4, 4 )

@test sum(s4_4_gt1 - s4_4.certainidx) == 0
@test length(s4_4.allelements[1]) > 0.5*s4_4_gt2[2]
@test sum(s4_4_gt2[2] .- s4_4.allelements[2]) == 0
@test length(s4_4.allelements[3]) > s4_4_gt2[2]
@test length(s4_4.allelements[4]) > s4_4_gt2[3]
@test length(s4_4.allelements[5]) > s4_4_gt2[4]
@test length(s4_4.allelements[1]) + length(s4_4.allelements[3]) + length(s4_4.allelements[4]) + length(s4_4.allelements[5]) == s4_4_gt2[5]

@test s4_4_gt3[1][1] == s4_4.activehypo[1][1]
@test s4_4_gt3[2][1] == s4_4.activehypo[2][1]
@test s4_4_gt3[3][1] == s4_4.activehypo[3][1]
@test s4_4_gt3[4][1] == s4_4.activehypo[4][1]
@test sum(s4_4_gt3[1][2] .- s4_4.activehypo[1][2]) == 0
@test sum(s4_4_gt3[2][2] .- s4_4.activehypo[2][2]) == 0
@test sum(s4_4_gt3[3][2] .- s4_4.activehypo[3][2]) == 0
@test sum(s4_4_gt3[4][2] .- s4_4.activehypo[4][2]) == 0

@test sum(s4_4.mhidx .== 0) > s4_4_gt2[2]
@test sum(s4_4.mhidx .== 2) > s4_4_gt2[2]
@test sum(s4_4.mhidx .== 3) > s4_4_gt2[3]
@test sum(s4_4.mhidx .== 4) > s4_4_gt2[4]

@test sum( [1:N;][s4_4.mhidx .== 0] .== s4_4.allelements[1] ) == length(s4_4.allelements[1])
@test 0 == length(s4_4.allelements[2])
@test sum( [1:N;][s4_4.mhidx .== 2] .== s4_4.allelements[3] ) == length(s4_4.allelements[3])
@test sum( [1:N;][s4_4.mhidx .== 3] .== s4_4.allelements[4] ) == length(s4_4.allelements[4])
@test sum( [1:N;][s4_4.mhidx .== 4] .== s4_4.allelements[5] ) == length(s4_4.allelements[5])
@test length(s4_4.mhidx) == s4_4_gt4

@warn "only partially testing tri-modality"

##

end
