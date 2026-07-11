# Simple examples to illustrate how to obtain a Bayes (junction) tree using
# with beautiful Tex labels in `.dot` and `.tex` format.

using Test

using IncrementalInference

# EXPERIMENTAL FEATURE, 4Q19: need `sudo apt install dot2tex`
import IncrementalInference: generateTexTree

@testset "testing generateTexTree" begin

# Create a dummy factor graph, with variables and constraints.
fg = initfg()

# Add four pose variables, with 'x' symbol.
addVariable!(fg, :x1, ContinuousScalar)
addVariable!(fg, :x2, ContinuousScalar)
addVariable!(fg, :x3, ContinuousScalar)
addVariable!(fg, :x4, ContinuousScalar)

# Add two landmark variables, with 'l' symbol.
addVariable!(fg, :l1, ContinuousScalar)
addVariable!(fg, :l2, ContinuousScalar)

# Add the pose chain constraints (odometry and priors).
addFactor!(fg, [:x1], Prior(Normal()))
addFactor!(fg, [:x1;:x2], LinearRelative(Normal()))
addFactor!(fg, [:x2;:x3], LinearRelative(Normal()))
addFactor!(fg, [:x3;:x4], LinearRelative(Normal()))

# Add the pose-landmark constraints (range measurements)
addFactor!(fg, [:x1;:l1], LinearRelative(Normal()))
addFactor!(fg, [:x2;:l1], LinearRelative(Normal()))
addFactor!(fg, [:x3;:l1], LinearRelative(Normal()))
addFactor!(fg, [:x2;:l2], LinearRelative(Normal()))
addFactor!(fg, [:x3;:l2], LinearRelative(Normal()))
addFactor!(fg, [:x4;:l2], LinearRelative(Normal()))

# Let's take a peek to see how our factor graph looks like.
# drawGraph(fg, show=true)
# As well as our tree (AMD ordering)
tree = buildTreeReset!(fg)
# drawTree(tree, show=true, imgs=false)

# Now, let's generate the corresponding `.dot` and `.tex`.
texTree = generateTexTree(tree, filepath=joinpath(@__DIR__,"tmp","bt"))

# All you have to do now is compile your newly created `.tex` file, probably
# include the `bm` package (`\usepackage{bm}`), and enjoy!

# fake a count of 1, since we are actually testing generateTexTree
@test true

end  # testset
