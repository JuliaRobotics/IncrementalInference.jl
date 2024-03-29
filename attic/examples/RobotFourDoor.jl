# Example file showing the classic 1D robot and four door example
# Best if run piece by piece in a Julia REPL or VSCode IDE


using IncrementalInference

# [OPTIONAL] plotting libraries are loaded separately
using Cairo, RoMEPlotting
Gadfly.set_default_plot_size(35cm,20cm)


## example parameters
# Number of kernels representing each marginal belief
N=100

# prior knowledge of four possible door locations
cv = 3.0
doorPrior = Mixture(Prior,
                    [Normal(-100,cv);Normal(0,cv);Normal(100,cv);Normal(300,cv)],
                    [1/4;1/4;1/4;1/4] )

## Build the factor graph object
fg = initfg()
 # WIP, will become default option in the future
getSolverParams(fg).useMsgLikelihoods = true

# first pose location
v1 = addVariable!(fg,:x1,ContinuousScalar,N=N)
# see a door for the first time
addFactor!(fg,[:x1], doorPrior)

# first solution with only one variable and factor (may take a few moments on first JIT compiling)
solveTree!(fg)
plotKDE(fg, :x1)


## drive to second pose location
addVariable!(fg,:x2, ContinuousScalar, N=N)
addFactor!(fg,[:x1;:x2],LinearRelative(Normal(50.0,2.0)))

# drive to third pose location
v3=addVariable!(fg,:x3,ContinuousScalar, N=N)
addFactor!(fg,[:x2;:x3], LinearRelative( Normal(50.0,4.0)))

# see a door for the second time
addFactor!(fg,[:x3], doorPrior)

# second solution should be much quicker
solveTree!(fg)
plotKDE(fg, [:x1; :x2; :x3])

# drive to forth and final pose location
addVariable!(fg,:x4,ContinuousScalar, N=N)
addFactor!(fg,[:x3;:x4], LinearRelative( Normal(200.0,4.0)))

# lets see the prediction of where pose :x4 might be
initAll!(fg)
plotKDE(fg, :x4)


## make a third door sighting
addFactor!(fg,[:x4], doorPrior)

# solve over all data
tree = solveTree!(fg)

# list variables and factors in fg
@show ls(fg) # |> sortDFG
@show lsf(fg)

## draw all beliefs
pl = plotKDE(fg, [:x1;:x2;:x3;:x4])

#save plot to file
pl |> PNG("4doors.png") # can also do SVG, PDF


## If interested, here is the Bayes/Junction tree too
# drawTree(tree, show=true) # using Graphviz and Linux evince for pdf


#
