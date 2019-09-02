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

writeGraphPdf(fg, show=true)

eo = [:x1; :l3; :l1; :x5; :x2; :l2; :x4; :x3]

tree = buildTreeFromOrdering!(fg,eo)

# using Gadfly, Fontconfig, Cairo
# drawTree(tree, show=true, imgs=false)
# drawTree(tree, show=true, imgs=true)
#


end



## Comparible reference available from gtsam/testSymbolicEliminationTree.cpp

#        l1                  l2
#    /    |                /  |
# x1 --- x2 --- x3 --- x4 --- x5
#                          \  |
#                            l3

# (X(1), L(1));
# (X(1), X(2));
# (X(2), L(1));
# (X(2), X(3));
# (X(3), X(4));
# (X(4), L(2));
# (X(4), X(5));
# (L(2), X(5));
# (X(4), L(3));
# (X(5), L(3));

# EliminationTreeTester::MakeTree(list_of
#     (MakeNode(X(3), SymbolicFactorGraph(), list_of
#       (MakeNode(X(2), list_of(SymbolicFactor(X(2), X(3))), list_of
#         (MakeNode(L(1), list_of(SymbolicFactor(X(2), L(1))), list_of
#           (MakeNode(X(1), list_of(SymbolicFactor(X(1), L(1))) (SymbolicFactor(X(1), X(2)))))))))
#       (MakeNode(X(4), list_of(SymbolicFactor(X(3), X(4))), list_of
#         (MakeNode(L(2), list_of(SymbolicFactor(X(4), L(2))), list_of
#           (MakeNode(X(5), list_of(SymbolicFactor(X(4), X(5))) (SymbolicFactor(L(2), X(5))), list_of
#             (MakeNode(L(3), list_of(SymbolicFactor(X(4), L(3))) (SymbolicFactor(X(5), L(3))))))))))) )));


#
