# Example file showing the classic 1D robot and four door example
# Best if run piece by piece in a Julia REPL or VSCode IDE


using IncrementalInference

# [OPTIONAL] plotting libraries are loaded separately
using Cairo, RoMEPlotting
Gadfly.set_default_plot_size(35cm,20cm)


## Ground truth and parameters
gt = Dict{Symbol, Matrix{Float64}}()
# HMM computed ground truth for first 3 poses only
gt[:x1]=[-100.0 0.0; 1.96 1.96]
gt[:x2]=[-50.0 50; 3.1 3.1]
gt[:x3]=[100.0 0.0; 3.05 3.05]

# Number of kernels representing each marginal belief
N=100

# prior knowledge of four possible door locations
cv = 3.0
doorPrior = Mixture(Prior,
                    [Normal(-100,cv);Normal(0,cv);Normal(100,cv);Normal(300,cv)],
                    [1/4;1/4;1/4;1/4] )

## Build the factor graph object
fg = initfg()
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
ensureAllInitialized!(fg)
plotKDE(fg, :x4)


## make a third door sighting
addFactor!(fg,[:x4], doorPrior)

# solve over all data
tree, smt, hists = solveTree!(fg)

# list variables and factors in fg
@show ls(fg)
@show lsf(fg)

## draw all beliefs
pl = plotKDE(fg, [:x1;:x2;:x3;:x4])

#save plot to file
pl |> PNG("4doors.png") # can also do SVG, PDF


## If interested, here is the Bayes/Junction tree too
# drawTree(tree, show=true) # using Graphviz and Linux evince for pdf











# if false
#
#   # draw upward messages
#   msgPlots=drawTreeUpwardMsgs(fg, tree, N=500); # to init memory for eval(parse(string))
#   println("Upward messages for all cliques except root")
#   # vvMsgs = vstackedDensities(fg, tree, msgPlots)
#   Gadfly.set_default_plot_size(17cm, 15cm)
#   Gadfly.draw(PNG("results/testMsgs.png",17cm,15cm),vvMsgs)
#   # vvMsgs
# end
