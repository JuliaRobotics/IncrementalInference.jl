using Test
using Distributed

addprocs(2)
using IncrementalInference
@everywhere using IncrementalInference

@testset "test multiprocess solveTree!" begin

fg = generateCanonicalFG_Kaess()

solveTree!(fg)

end
