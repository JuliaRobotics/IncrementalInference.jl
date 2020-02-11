# basic graph tests for trivial cases

using Test
using Statistics
using IncrementalInference


@testset "Test basic single variable graph with one prior..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(0.0,1.0)))

tree, smt, hist = solveTree!(fg)

# check mean and covariance
@test (getKDE(fg, :x0) |> getKDEMean .|> abs)[1] < 0.5
@test 0.4 < Statistics.cov( getPoints(getKDE(fg, :x0))[1,:] ) < 1.8

end


@testset "Test basic single variable graph with one prior offset by 1000..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(1000.0,1.0)))

tree, smt, hist = solveTree!(fg)

# check mean and covariance
@test abs((getKDE(fg, :x0) |> getKDEMean)[1]-1000) < 0.5
@test 0.4 < Statistics.cov( getPoints(getKDE(fg, :x0))[1,:] ) < 1.8

end


@testset "Test basic single variable graph with two identical priors..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(0.0,1.0)))
addFactor!(fg, [:x0;], Prior(Normal(0.0,1.0)))

tree, smt, hist = solveTree!(fg)

# check mean and covariance
@test (getKDE(fg, :x0) |> getKDEMean .|> abs)[1] < 0.3
# should be sqrt(0.5) = 0.7, but lands near 0.6 instead -- computation is too confident.
@test 0.3 < Statistics.cov( getPoints(getKDE(fg, :x0))[1,:] ) < 1.0

end


@testset "Test basic single variable graph with three identical priors..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(0.0,1.0)))
addFactor!(fg, [:x0;], Prior(Normal(0.0,1.0)))
addFactor!(fg, [:x0;], Prior(Normal(0.0,1.0)))

tree, smt, hist = solveTree!(fg)

# check mean and covariance
@test (getKDE(fg, :x0) |> getKDEMean .|> abs)[1] < 0.4
# should be sqrt(1/3) = 0.577, but lands near 0.35 instead -- computation is too confident.
@test 0.1 < Statistics.cov( getPoints(getKDE(fg, :x0))[1,:] ) < 0.75

end



@testset "Test basic single variable graph with two priors at + and - 1..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(-1.0,1.0)))
addFactor!(fg, [:x0;], Prior(Normal(+1.0,1.0)))

tree, smt, hist = solveTree!(fg)

# check mean and covariance -- should be zero
@test (getKDE(fg, :x0) |> getKDEMean .|> abs)[1] < 0.8
# should be sqrt(1/2) = 0.707 -- computation results nearer 0.7.
@test 0.2 < Statistics.cov( getPoints(getKDE(fg, :x0))[1,:] ) < 1

end


@testset "Test basic single variable graph with two priors at + and - 1, offset by -1000..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(-1.0-1000,1.0)))
addFactor!(fg, [:x0;], Prior(Normal(+1.0-1000,1.0)))

tree, smt, hist = solveTree!(fg)

# check mean and covariance -- should be zero
@test abs((getKDE(fg, :x0) |> getKDEMean)[1] + 1000) < 0.6
# should be sqrt(1/2) = 0.707 -- computation results nearer 0.7.
@test 0.2 < Statistics.cov( getPoints(getKDE(fg, :x0))[1,:] ) < 1

end



@testset "Test basic two variable graph with two identical priors and weak connection..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(0.0,1.0)))
addFactor!(fg, [:x1;], Prior(Normal(0.0,1.0)))
addFactor!(fg, [:x0;:x1;], LinearConditional(Normal(0.0,10.0)))

tree, smt, hist = solveTree!(fg)

# check mean and covariance -- should be zero
@test (getKDE(fg, :x0) |> getKDEMean .|> abs)[1] < 0.6
@test (getKDE(fg, :x1) |> getKDEMean .|> abs)[1] < 0.6

@test 0.4 < Statistics.cov( getPoints(getKDE(fg, :x0))[1,:] ) < 2.2
@test 0.4 < Statistics.cov( getPoints(getKDE(fg, :x1))[1,:] ) < 2.3

end


@testset "Test basic two variable graph with two separated priors and weak connection..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(-1.0,1.0)))
addFactor!(fg, [:x1;], Prior(Normal(+1.0,1.0)))
addFactor!(fg, [:x0;:x1;], LinearConditional(Normal(0.0,10.0)))

tree, smt, hist = solveTree!(fg)

# check mean and covariance -- should be near each prior
@test abs((getKDE(fg, :x0) |> getKDEMean)[1]+1) < 0.5
@test abs((getKDE(fg, :x1) |> getKDEMean)[1]-1) < 0.5

@test 0.4 < Statistics.cov( getPoints(getKDE(fg, :x0))[1,:] ) < 2.5
@test 0.4 < Statistics.cov( getPoints(getKDE(fg, :x1))[1,:] ) < 2.5

end


@testset "Test basic two variable graph with two separated priors and strong connection..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)
addVariable!(fg, :x2, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal(-1.0,1.0)))
addFactor!(fg, [:x2;], Prior(Normal(+1.0,1.0)))
addFactor!(fg, [:x0;:x1;], LinearConditional(Normal(0.0,1.0)))
addFactor!(fg, [:x1;:x2;], LinearConditional(Normal(0.0,1.0)))

tree, smt, hist = solveTree!(fg)



# check mean and covariance -- should between two priors somewhere
@test abs((getKDE(fg, :x0) |> getKDEMean)[1] + 1) < 0.8
@test abs((getKDE(fg, :x1) |> getKDEMean)[1]) < 0.9
@test abs((getKDE(fg, :x2) |> getKDEMean)[1] - 1) < 0.9

@test 0.3 < Statistics.cov( getPoints(getKDE(fg, :x0))[1,:] ) < 1.8
@test 0.3 < Statistics.cov( getPoints(getKDE(fg, :x1))[1,:] ) < 1.8
@test 0.3 < Statistics.cov( getPoints(getKDE(fg, :x2))[1,:] ) < 1.8

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
addFactor!(fg, [:x0;:x1;], LinearConditional(Normal(0.0,1.0)))
addFactor!(fg, [:x1;:x2;], LinearConditional(Normal(0.0,1.0)))
addFactor!(fg, [:x2;:x3;], LinearConditional(Normal(0.0,1.0)))
addFactor!(fg, [:x3;:x4;], LinearConditional(Normal(0.0,1.0)))

# drawGraph(fg, show=true)

tree, smt, hist = solveTree!(fg)

# using KernelDensityEstimatePlotting
# plotKDE((x->getKDE(fg,x)).([:x0;:x1;:x2;:x3;:x4]))
# using Gadfly, Cairo, Fontconfig
# drawTree(tree,show=true,imgs=true)

# check mean and covariance -- should be zero
X0 = (getKDE(fg, :x0) |> getKDEMean)[1]
X1 = (getKDE(fg, :x1) |> getKDEMean)[1]
X2 = (getKDE(fg, :x2) |> getKDEMean)[1]
X3 = (getKDE(fg, :x3) |> getKDEMean)[1]
X4 = (getKDE(fg, :x4) |> getKDEMean)[1]

@test X0 < X1 < X2 < X3 < X4

@test abs(X0+X4) < 2.2
@test abs(X1+X3) < 2.2
@test abs(X2) < 2.2

@test 0.2 < Statistics.cov( getPoints(getKDE(fg, :x0))[1,:] ) < 2.3
@test 0.2 < Statistics.cov( getPoints(getKDE(fg, :x1))[1,:] ) < 2.3
@test 0.2 < Statistics.cov( getPoints(getKDE(fg, :x2))[1,:] ) < 2.4
@test 0.2 < Statistics.cov( getPoints(getKDE(fg, :x3))[1,:] ) < 2.4
@test 0.2 < Statistics.cov( getPoints(getKDE(fg, :x4))[1,:] ) < 2.5

end










#
