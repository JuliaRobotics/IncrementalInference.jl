using Test
using NLsolve
using Distributions
# using Statistics

using IncrementalInference

import IncrementalInference: getSample


##

@testset "Ensure lambda's work with NLsolve" begin

function f!(F, x, dummy)
    F[1] = (x[1]+3)*(x[2]^3-7)+18
    F[2] = sin(x[2]*exp(x[1])-1)
end

data = randn(2,3)
r = nlsolve( (res, x) -> f!(res, x, data), zeros(2), inplace=true )

@test norm(r.zero - [0, 1.0]) < 1e-10

end

#
