# lets test long approxConv chains


using Test
using IncrementalInference
using Statistics
using TensorCast

##

@testset "zdim size of factor is correct" begin

##

fg = initfg()
addVariable!(fg, :x0, ContinuousEuclid{3})
addVariable!(fg, :x1, ContinuousEuclid{3})

lr3 =  LinearRelative(MvNormal([1;0;0.0],diagm(ones(3))))
addFactor!(fg, [:x0; :x1], lr3, graphinit=false)

ccw = IIF._getCCW(fg, :x0x1f1)

@test IIF._getZDim(ccw) == 3

res = testFactorResidualBinary( lr3, ContinuousEuclid{3}, ContinuousEuclid{3},
                                zeros(3),[0;0;1.0], ([0;0;0.5],))
#

@test (sum(abs.(res)) - 0.5) < 1e-10


##

end


@testset "approxConv basic and chains" begin

##

fg = generateCanonicalFG_Kaess()

# from a prior to neighbor
pts_ = approxConv(fg, :x1f1, :x1)

# lazy
@cast pts[i,j] := pts_[j][i]

@test Statistics.mean(pts) |> abs  < 0.4
@test 0.5 < Statistics.std(pts) < 1.5

# set a value in graph to start things off
initManual!(fg, :x1, pts_)

# legacy case where relative to neighbor 
pts_ = approxConv(fg, :x1x2f1, :x2)
@cast pts[i,j] := pts_[j][i]


@test Statistics.mean(pts) |> abs  < 0.7
@test 0.7 < Statistics.std(pts) < 2

# along a chain of variables
pts_ = approxConv(fg, :x1, :x3)
@cast pts[i,j] := pts_[j][i]

@test Statistics.mean(pts) |> abs  < 1.5
@test 1.3 < Statistics.std(pts) < 3.0

# from a prior down the chain of variables
pts_ = approxConv(fg, :x1f1, :l2)
@cast pts[i,j] := pts_[j][i]

@test Statistics.mean(pts) |> abs  < 1.5
@test 1.6 < Statistics.std(pts) < 4.0

##

end



#