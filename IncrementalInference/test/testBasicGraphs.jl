# basic graph tests for trivial cases

using Test
using Statistics
using IncrementalInference

##

@testset "test basic utility functions" begin

@test incrSuffix(:x45_4) == :x45_5

@test incrSuffix(:x45, +3) == :x48

@test incrSuffix(:x45_4, -1) == :x45_3

end

@testset "Test basic single variable graph with one prior..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(0.0,1.0)))

# test solved flag
@test getSolvedCount(fg, :x0) == 0
@test !isSolved(getVariable(fg, :x0))

makeSolverData!(fg)

# run solver once
tree = solveTree!(fg)

@test getSolvedCount(fg, :x0) == 1
@test isSolved(fg, :x0)

tree = solveTree!(fg)

@test getSolvedCount(fg, :x0) == 2
@test isSolved(fg, :x0)

# check mean and covariance
@test (getBelief(fg, :x0) |> getKDEMean .|> abs)[1] < 0.5
pts_ = getPoints(getBelief(fg, :x0))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.3 < Statistics.cov( pts[1,:] ) < 1.9

# test free solvable variables (occurs in fixed-/ clique recycling)
addVariable!(fg, :x1, ContinuousScalar, solvable=1)

solveTree!(fg, storeOld=true)

@test getSolvable(fg, :x1) == 0

end


@testset "Test basic single variable graph with one prior offset by 1000..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(1000.0,1.0)))

tree = solveTree!(fg)

# check mean and covariance
@test abs((getBelief(fg, :x0) |> getKDEMean)[1]-1000) < 0.5
pts_ = getPoints(getBelief(fg, :x0))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.4 < Statistics.cov( pts[1,:] ) < 1.8

end


@testset "Test basic single variable graph with two identical priors..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(0.0,1.0)))
addFactor!(fg, [:x0;], Prior(Normal(0.0,1.0)))

tree = solveTree!(fg)

# check mean and covariance
@test (getBelief(fg, :x0) |> getKDEMean .|> abs)[1] < 0.4
# should be sqrt(0.5) = 0.7, but lands near 0.6 instead -- computation is too confident.
pts_ = getPoints(getBelief(fg, :x0))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.3 < Statistics.cov( pts[1,:] ) < 1.0

end


@testset "Test basic single variable graph with three identical priors..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(0.0,1.0)))
addFactor!(fg, [:x0;], Prior(Normal(0.0,1.0)))
addFactor!(fg, [:x0;], Prior(Normal(0.0,1.0)))

tree = solveTree!(fg)

# check mean and covariance
@test (getBelief(fg, :x0) |> getKDEMean .|> abs)[1] < 0.4
# should be sqrt(1/3) = 0.577, but lands near 0.35 instead -- computation is too confident.
pts_ = getPoints(getBelief(fg, :x0))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.1 < Statistics.cov( pts[1,:] ) < 0.75

end



@testset "Test basic single variable graph with two priors at + and - 1..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(-1.0,1.0)))
addFactor!(fg, [:x0;], Prior(Normal(+1.0,1.0)))

tree = solveTree!(fg)

# check mean and covariance -- should be zero
@test (getBelief(fg, :x0) |> getKDEMean .|> abs)[1] < 0.8
# should be sqrt(1/2) = 0.707 -- computation results nearer 0.7.
pts_ = getPoints(getBelief(fg, :x0))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.2 < Statistics.cov( pts[1,:] ) < 1.5

end


@testset "Test basic single variable graph with two priors at + and - 1, offset by -1000..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(-1.0-1000,1.0)))
addFactor!(fg, [:x0;], Prior(Normal(+1.0-1000,1.0)))

tree = solveTree!(fg)

# check mean and covariance -- should be zero
@test abs((getBelief(fg, :x0) |> getKDEMean)[1] + 1000) < 0.6
# should be sqrt(1/2) = 0.707 -- computation results nearer 0.7.
pts_ = getPoints(getBelief(fg, :x0))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.2 < Statistics.cov( pts[1,:] ) < 1.1

end



@testset "Test basic two variable graph with two identical priors and weak connection..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(0.0,1.0)))
addFactor!(fg, [:x1;], Prior(Normal(0.0,1.0)))
addFactor!(fg, [:x0;:x1;], LinearRelative(Normal(0.0,10.0)))

tree = solveTree!(fg)

# check mean and covariance -- should be zero
@test (getBelief(fg, :x0) |> getKDEMean .|> abs)[1] < 0.6
@test (getBelief(fg, :x1) |> getKDEMean .|> abs)[1] < 0.6

pts_ = getPoints(getBelief(fg, :x0))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.4 < Statistics.cov( pts[1,:] ) < 2.3
pts_ = getPoints(getBelief(fg, :x1))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.4 < Statistics.cov( pts[1,:] ) < 2.4

end


@testset "Test basic two variable graph with two separated priors and weak connection..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(-1.0,1.0)))
addFactor!(fg, [:x1;], Prior(Normal(+1.0,1.0)))
addFactor!(fg, [:x0;:x1;], LinearRelative(Normal(0.0,10.0)))

tree = solveTree!(fg)

# check mean and covariance -- should be near each prior
@test abs((getBelief(fg, :x0) |> getKDEMean)[1]+1) < 0.75
@test abs((getBelief(fg, :x1) |> getKDEMean)[1]-1) < 0.75

pts_ = getPoints(getBelief(fg, :x0))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.3 < Statistics.cov( pts[1,:] ) < 2.5

pts_ = getPoints(getBelief(fg, :x1))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.3 < Statistics.cov( pts[1,:] ) < 2.5

end


@testset "Test basic two variable graph with two separated priors and strong connection..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)
addVariable!(fg, :x2, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(-1.0,1.0)))
addFactor!(fg, [:x2;], Prior(Normal(+1.0,1.0)))
addFactor!(fg, [:x0;:x1;], LinearRelative(Normal(0.0,1.0)))
addFactor!(fg, [:x1;:x2;], LinearRelative(Normal(0.0,1.0)))

tree = solveTree!(fg)



# check mean and covariance -- should between two priors somewhere
@test abs((getBelief(fg, :x0) |> getKDEMean)[1] + 1) < 0.9
@test abs((getBelief(fg, :x1) |> getKDEMean)[1]) < 0.9
@test abs((getBelief(fg, :x2) |> getKDEMean)[1] - 1) < 0.9

pts_ = getPoints(getBelief(fg, :x0))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.3 < Statistics.cov( pts[1,:] ) < 1.8

pts_ = getPoints(getBelief(fg, :x1))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.3 < Statistics.cov( pts[1,:] ) < 2.0

pts_ = getPoints(getBelief(fg, :x2))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.3 < Statistics.cov( pts[1,:] ) < 2.2

end



@testset "Test basic five variable graph with two separated priors and nominal connection..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)
addVariable!(fg, :x2, ContinuousScalar)
addVariable!(fg, :x3, ContinuousScalar)
addVariable!(fg, :x4, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(-3.0,1.0)))
addFactor!(fg, [:x4;], Prior(Normal(+3.0,1.0)))
addFactor!(fg, [:x0;:x1;], LinearRelative(Normal(0.0,1.0)))
addFactor!(fg, [:x1;:x2;], LinearRelative(Normal(0.0,1.0)))
addFactor!(fg, [:x2;:x3;], LinearRelative(Normal(0.0,1.0)))
addFactor!(fg, [:x3;:x4;], LinearRelative(Normal(0.0,1.0)))

# #1196
drawGraph(fg, filepath="testgraphplot/myfg.dot", show=false)

tree = solveTree!(fg, storeOld=true)

# using KernelDensityEstimatePlotting
# plotKDE((x->getBelief(fg,x)).([:x0;:x1;:x2;:x3;:x4]))
# using Gadfly, Cairo, Fontconfig
# drawTree(tree,show=true,imgs=true)

# check mean and covariance -- should be zero
X0 = (getBelief(fg, :x0) |> getKDEMean)[1]
X1 = (getBelief(fg, :x1) |> getKDEMean)[1]
X2 = (getBelief(fg, :x2) |> getKDEMean)[1]
X3 = (getBelief(fg, :x3) |> getKDEMean)[1]
X4 = (getBelief(fg, :x4) |> getKDEMean)[1]

@test X0 < X1 < X2 < X3 < X4

@test abs(X0+X4) < 2.2
@test abs(X1+X3) < 2.2
@test abs(X2) < 2.2

pts_ = getPoints(getBelief(fg, :x0))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.2 < Statistics.cov( pts[1,:] ) < 2.8

pts_ = getPoints(getBelief(fg, :x1))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.2 < Statistics.cov( pts[1,:] ) < 2.9

pts_ = getPoints(getBelief(fg, :x2))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.2 < Statistics.cov( pts[1,:] ) < 3.0

pts_ = getPoints(getBelief(fg, :x3))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.2 < Statistics.cov( pts[1,:] ) < 3.1

pts_ = getPoints(getBelief(fg, :x4))
TensorCast.@cast pts[i,j] := pts_[j][i]
@test 0.2 < Statistics.cov( pts[1,:] ) < 3.2



@testset "Test localProduct on solveKey" begin

localProduct(fg,:x2)
localProduct(fg,:x2, solveKey=:graphinit)

end

end


##
@testset "consistency check on more factors (origin is a DERelative fail case)" begin
##

fg = initfg()

addVariable!(fg, :x0, Position{1})
addFactor!(fg, [:x0], Prior(Normal(1.0, 0.01)))

# force a basic setup
initAll!(fg)
@test isapprox( 1, calcMeanMaxSuggested(fg, :x0, :default).suggested[1]; atol=0.1)

##

addVariable!(fg, :x1, Position{1})
addFactor!(fg, [:x0;:x1], LinearRelative(Normal(1.0, 0.01)))

addVariable!(fg, :x2, Position{1})
addFactor!(fg, [:x1;:x2], LinearRelative(Normal(1.0, 0.01)))

addVariable!(fg, :x3, Position{1})
addFactor!(fg, [:x2;:x3], LinearRelative(Normal(1.0, 0.01)))

##

tree = solveGraph!(fg)

##

@test isapprox( 1, calcMeanMaxSuggested(fg, :x0, :default).suggested[1]; atol=0.1)
@test isapprox( 4, calcMeanMaxSuggested(fg, :x3, :default).suggested[1]; atol=0.3)

## check contents of tree messages

tree[1]
msg1 = IIF.getMessageBuffer(tree[1])

##
end



@testset "Test graph reset to init..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(1000.0,1.0)))

initAll!(fg)

# init values before solve
pts_ = getPoints(getBelief(fg, :x0))
TensorCast.@cast pts[i,j] := pts_[j][i]
X0 = pts |> deepcopy

tree = solveTree!(fg)

# values after solve
pts_ = getPoints(getBelief(fg, :x0))
TensorCast.@cast pts[i,j] := pts_[j][i]
X0s = pts

@test 1e-10 < norm(X0 - X0s)

resetInitialValues!(fg)

pts_ = getPoints(getBelief(fg, :x0))
TensorCast.@cast pts[i,j] := pts_[j][i]
X0reset = pts |> deepcopy

@test norm(X0 - X0reset) < 1e-10


end


@testset "Test MetaPrior" begin

fg = generateGraph_Kaess()

addFactor!(fg, [:x1], MetaPrior(nothing))

solveGraph!(fg)
IIF.solveGraphParametric!(fg)

end


#
