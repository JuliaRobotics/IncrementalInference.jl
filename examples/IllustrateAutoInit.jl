# example to illustate how autoinit functions work, and used as a development script towards a standard unit test.

using IncrementalInference

import IncrementalInference: getSample

## MARKER START is a required block of code, but can be reviewed at the end of the tutorial,
# Run this but skip past these user defined functions for a quicker introduction======================


# and a bi-modal conditional function

## TODO: NEW STANDARD FEATURE USE addFactor!(fg, .., multihypo=[...]),   or  MixtureLinearConditional
# struct MixtureConditional <: IncrementalInference.FunctorPairwise
#   z::Vector{Distribution}
#   c::Categorical
# end
# getSample(s::MixtureConditional, N::Int=1) = (rand.(s.z, N)..., rand(s.c, N))
# function (s::MixtureConditional)(res::Array{Float64},
#       userdata::FactorMetadata,
#       idx::Int,
#       meas::Tuple,
#       X1::Array{Float64,2},
#       X2::Array{Float64,2}
#   )
#   res[1] = meas[meas[end][idx]][idx] - (X2[1,idx] - X1[1,idx])
#   nothing
# end

## Define a model with these user defined functions ===============================================
## MARKER END



# Start with an empty graph
fg = initfg()

# add the first node
addVariable!(fg, :x0, ContinuousScalar)

# this is unary (prior) factor and does not immediately trigger autoinit of :x0.
addFactor!(fg, [:x0], Prior(Normal(0,1)))

# writeGraphPdf(fg, file="/home/dehann/Downloads/fgx0.png")

# To visualize the factor graph structure, assuming GraphViz.jl is available.
# Graphs.plot(fg.g)
# Please find `writeGraphPdf` member definition at the end of this tutorial
# writeGraphPdf(fg, file="/home/dehann/Downloads/fgx01.png")

# Also consider automatic initialization of variables.
@show isInitialized(fg, :x0)
# Why is x0 not initialized?
# Since no other variable nodes have been 'connected to' (or depend) on :x0,
# and future intentions of the user are unknown, the initialization of :x0 is
# deferred until such criteria are met.

# Auto initialization of :x0 is triggered by connecting the next (probably uninitialized) variable node
addVariable!(fg, :x1, ContinuousScalar)

# with a linear conditional belief to :x0
# P(Z | :x1 - :x0 ) where Z ~ Normal(10,1)
addFactor!(fg, [:x0, :x1], LinearConditional(Normal(10.0,1)))

# writeGraphPdf(fg, file="/home/dehann/Downloads/fgx01.png")

# x0 should not be initialized
@show isInitialized(fg, :x0)

# To plot the belief states of variables -- (assuming `RoMEPlotting` is available)
# Remember the first time executions are slow given required code compilation,
# and that future versions of these package will use more precompilation
# to reduce first execution running cost.
using RoMEPlotting

plotKDE(fg, :x0)

# Since no other variables 'are yet connected'/depend on :x1, it will not be initialized
@show isInitialized(fg, :x1)

# we can force all the variable nodes to initialize
ensureAllInitialized!(fg)

# now draw both :x0
plotKDE(fg, [:x0, :x1])

# add another node, but introduce more general beliefs
addVariable!(fg, :x2, ContinuousScalar)

mmo = MixtureConditional([Rayleigh(3); Uniform(30,55)], Categorical([0.4; 0.6]))
addFactor!(fg, [:x1, :x2], mmo)

# Graphs.plot(fg.g)
# writeGraphPdf(fg, file="/home/dehann/Downloads/fgx012.png")

# By again forcing the initialization of :x3 for illustration
ensureAllInitialized!(fg)

# the predicted marginal probability densities are
plotKDE(fg, [:x0, :x1, :x2])

# Now transmit this 'weird' multi-modal marginal belief through a another unimodal linear offset (conditional likelihood)
addVariable!(fg, :x3, ContinuousScalar)

addFactor!(fg, [:x2, :x3], LinearConditional(Normal(-50, 1)))
# note, this addFactor step relies on :x2 being initialized and would have done so if we didn't call ensureAllInitialized! a few lines earlier.

# writeGraphPdf(fg, file="/home/dehann/Downloads/fgx0123.png")

ensureAllInitialized!(fg)
plotKDE(fg, [:x0, :x1, :x2, :x3])


lo3 = LinearConditional(Normal(40, 1))
addFactor!(fg, [:x3, :x0], lo3)


# writeGraphPdf(fg, file="/home/dehann/Downloads/fgx0123c.png")

plotKDE(fg, [:x0, :x1, :x2, :x3])



# Find global best likelihood solution (posterior belief)
# After defining the problem, we can find the 'minimum free energy' solution
tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)


# and look at the posterior belief, and notice which consensus modes stand out in the posterior
plotKDE(fg, [:x0, :x1, :x2, :x3])




# Helper function for graphing the factor graph structure (using GraphViz)


# using Gadfly
# pl.guides[1] = Gadfly.Guide.xlabel("")
# push!(pl.guides, Gadfly.Guide.ylabel("density"))
#
# Gadfly.draw(PNG("/home/dehann/Downloads/plx012.png", 10cm, 7cm), pl)

## should complete and add to RoMEPlotting
# import KernelDensityEstimatePlotting: plot
# import Gadfly: plot
# import Graphs: plot
# import RoMEPlotting: plot
# function plot(fgl::FactorGraph, sym::Symbol; api::DataLayerAPI=IncrementalInference.dlapi)
#   PX = getKDE(getVariable(fgl, sym))
#   plot(PX)
# end

#
