# Showcasing the available analysis tools for the Bayes (Junction) tree.
# using Revise

using IncrementalInference
using DistributedFactorGraphs # For `isSolvable` function.
using Combinatorics # For creating the variable ordering `permutations`.
using SuiteSparse.CHOLMOD: SuiteSparse_long # For CCOLAMD constraints.
using Gadfly # For histogram and scatter plots.
Gadfly.set_default_plot_size(35cm,25cm)


latex_fonts = Theme(major_label_font="CMU Serif", major_label_font_size=16pt,
                    minor_label_font="CMU Serif", minor_label_font_size=14pt,
                    key_title_font="CMU Serif", key_title_font_size=12pt,
                    key_label_font="CMU Serif", key_label_font_size=10pt)
Gadfly.push_theme(latex_fonts)

# Get tree for each variable ordering in a factor graph.
fg = generateCanonicalFG_Kaess(graphinit=false)
all_trees = getAllTrees(deepcopy(fg))

# scores stores: (tree key ID, nnz, cost fxn 1, cost fxn 2).
unsorted_scores = Vector{Tuple{Int, Float64, Float64, Float64}}()
for key in keys(all_trees)
    e = all_trees[key] # (Bayes tree, var order, nnz
    tree = e[1] # Get the Bayes tree.
    cost1 = getTreeCost_01(tree)
    cost2 = getTreeCost_02(tree)
    push!(unsorted_scores, (key, e[3], cost1, cost2))
end

# Sort them to make sure the keys are in order.
scores = sort(unsorted_scores)

# Separate scores into vectors for plotting.
all_nnzs = (x->(x[2])).(scores)
costs_01 = (x->(x[3])).(scores)
costs_02 = (x->(x[4])).(scores)

min_ids_02 = findall(x->x == minimum(costs_02), costs_02)
max_ids_02 = findall(x->x == maximum(costs_02), costs_02)

min_ids_nnz = findall(x->x == minimum(all_nnzs), all_nnzs)
max_ids_nnz = findall(x->x == maximum(all_nnzs), all_nnzs)

# Find the intersection between best on both rubrics (lower left quadrant).
best_ids = findall(x->x in min_ids_02, min_ids_nnz)
# Find good factorizations but bad trees (upper left quadrant).
bad_trees_good_mats_ids = findall(x->x in max_ids_02, min_ids_nnz)
# Find good trees with bad matrix factorizations (lower right quadrant).
good_trees_bad_mats_ids = min_ids_02[findall(x->x == maximum(all_nnzs[min_ids_02]), all_nnzs[min_ids_02])]

# Get AMDs variable ordering.
amd_ordering = getEliminationOrder(fg)
amd_tree = buildTreeFromOrdering!(deepcopy(fg), amd_ordering)
amd_tree_nnz = nnzTree(amd_tree)
amd_tree_cost02 = getTreeCost_02(amd_tree)

# Get CCOLAMD variable ordering. First bring in CCOLAMD.
include(normpath(Base.find_package("IncrementalInference"), "..", "ccolamd.jl"))
A, varsym, fctsym = getAdjacencyMatrixSparse(fg)
colamd_ordering = varsym[Ccolamd.ccolamd(A)]
colamd_tree = buildTreeFromOrdering!(deepcopy(fg), colamd_ordering)
colamd_tree_nnz = nnzTree(colamd_tree)
colamd_tree_cost02 = getTreeCost_02(colamd_tree)

# Now add the iSAM2 constraint.
cons = zeros(SuiteSparse_long, length(A.colptr) - 1)
cons[findall(x->x == :x3, varsym)[1]] = 1 # NOTE(tonioteran) hardcoded for Kaess' example.
ccolamd_ordering = varsym[Ccolamd.ccolamd(A, cons)]
ccolamd_tree = buildTreeFromOrdering!(deepcopy(fg), ccolamd_ordering)
ccolamd_tree_nnz = nnzTree(ccolamd_tree)
ccolamd_tree_cost02 = getTreeCost_02(ccolamd_tree)

# Plot data points and underlying histogram.
bincnt = 20
layers = []
push!(layers, Gadfly.layer(x=[amd_tree_nnz],
                           y=[amd_tree_cost02],
                           Theme(default_color=colorant"green")))
push!(layers, Gadfly.layer(x=[colamd_tree_nnz],
                           y=[colamd_tree_cost02],
                           Theme(default_color=colorant"blue")))
push!(layers, Gadfly.layer(x=[ccolamd_tree_nnz],
                           y=[ccolamd_tree_cost02],
                           Theme(default_color=colorant"red")))
push!(layers, Gadfly.layer(x=all_nnzs,
                           y=costs_02,
                           Geom.hexbin(xbincount=bincnt, ybincount=bincnt)))
pl = Gadfly.plot(layers...,
            Guide.xlabel("Number of non zeros [int]"),
            Guide.ylabel("Tree cost [cfxn2]"),
            Guide.manual_color_key("",
             ["AMD", "COLAMD", "iSAM2"],
             ["green", "blue", "red"]))

img = SVG("vo_cost_canon_kaess.svg", 6inch, 6inch)
Gadfly.draw(img, pl)
