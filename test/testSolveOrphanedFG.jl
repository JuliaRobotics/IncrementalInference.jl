# test forest of graphs can solve with CSM,

using IncrementalInference
using Test



@testset "Test forest of orphaned graphs"

fg = initfg()
addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0;], Prior(Normal()))
addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0;:x1], LinearConditional(Normal()))

addVariable!(fg, :x10, ContinuousScalar)
addFactor!(fg, [:x10;], Prior(Normal()))
addVariable!(fg, :x11, ContinuousScalar)
addFactor!(fg, [:x10;:x11], LinearConditional(Normal()))

# dfgplot(fg)
# solve factor graph with two orphaned components
tree, smt, hist = solveTree!(fg)

# test tree will have two different root nodes
# ...

end
