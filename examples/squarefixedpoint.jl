# sqrt example

using Distributions
using KernelDensityEstimate, KernelDensityEstimatePlotting
using IncrementalInference
using Gadfly, DataFrames

# requried for overloading with new factors
import IncrementalInference: getSample


include(joinpath(@__DIR__, "SquareRootTypes.jl"))

## Direct example

fg = initfg()

addVariable!(fg, :x, ContinuousScalar)
# addVariable!(fg, :y, x0,N=N)

# TODO perhaps make Mixture instead for multiple computation
pts = rand(Distributions.Normal(4.0,0.05),100)   #;rand(Distributions.Normal(144.0,0.05),N)]
md = kde!(pts)

npx = NumbersPrior(md)
addVariable!(fg, :xy, ContinuousScalar)

addFactor!(fg, [:xy;], npx)

xty = Square(Distributions.Normal(0.0,0.01))
addFactor!(fg, [:x,:xy], xty)

# drawGraph(fg)

# initialize from prior
doautoinit!(fg, :xy)
# initialize any random numbers for the square root initial value
manualinit!(fg, :x, randn(1,100))

# find solution
tree, smt, hist = solveTree!(fg)

## plot the result
plotKDE(map(x->getKDE(fg,x), [:x; :xy]))
