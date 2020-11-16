# lets test long approxConv chains


using Test
using IncrementalInference
using Statistics


@testset "approxConv basic and chains" begin


fg = generateCanonicalFG_Kaess()

# from a prior to neighbor
pts = approxConv(fg, :x1f1, :x1)

@test Statistics.mean(pts) |> abs  < 0.3
@test 0.5 < Statistics.std(pts) < 1.5

# set a value in graph to start things off
initManual!(fg, :x1, pts)

# legacy case where relative to neighbor 
pts = approxConv(fg, :x1x2f1, :x2)

@test Statistics.mean(pts) |> abs  < 0.4
@test 0.7 < Statistics.std(pts) < 2

# along a chain of variables
pts = approxConv(fg, :x1, :x3)

@test Statistics.mean(pts) |> abs  < 0.5
@test 1.3 < Statistics.std(pts) < 2.5

# from a prior down the chain of variables
pts = approxConv(fg, :x1f1, :l2)

@test Statistics.mean(pts) |> abs  < 0.65
@test 1.7 < Statistics.std(pts) < 3.0


end



#