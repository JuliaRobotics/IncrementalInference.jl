# another tree test

using IncrementalInference
using Test

@testset "second tree formation test" begin


fg = initfg()

addVariable!(fg, :x1, ContinuousScalar)
addVariable!(fg, :x2, ContinuousScalar)
addVariable!(fg, :x3, ContinuousScalar)
addVariable!(fg, :x4, ContinuousScalar)
addVariable!(fg, :x5, ContinuousScalar)
addVariable!(fg, :l1, ContinuousScalar)
addVariable!(fg, :l2, ContinuousScalar)
addVariable!(fg, :l3, ContinuousScalar)

addFactor!(fg, [:x1;:l1], LinearConditional(Normal()))
addFactor!(fg, [:x1;:x2], LinearConditional(Normal()))
addFactor!(fg, [:x2;:l1], LinearConditional(Normal()))
addFactor!(fg, [:x2;:x3], LinearConditional(Normal()))
addFactor!(fg, [:x3;:x4], LinearConditional(Normal()))
addFactor!(fg, [:x4;:l2], LinearConditional(Normal()))
addFactor!(fg, [:x4;:x5], LinearConditional(Normal()))
addFactor!(fg, [:l2;:x5], LinearConditional(Normal()))
addFactor!(fg, [:x4;:l3], LinearConditional(Normal()))
addFactor!(fg, [:x5;:l3], LinearConditional(Normal()))

#writeGraphPdf(fg, show=true)

eo = [:x1; :l3; :l1; :x5; :x2; :l2; :x4; :x3]

tree = buildTreeFromOrdering!(fg,eo)

# using Gadfly, Fontconfig, Cairo
# drawTree(tree, show=true, imgs=false)
# drawTree(tree, show=true, imgs=true)
#
end


@testset "Variable ordering Bayes tree member check." begin
    # Would be nice to just be able to do that.
    fg = loadCanonicalFG_Kaess()
    # Choose specific variable ordering and perform check.
    vo = [:l1, :l2, :x1, :x2, :x3]
    tree = buildTreeFromOrdering!(fg, vo)
    @test vo == tree.variableOrder
end



#
