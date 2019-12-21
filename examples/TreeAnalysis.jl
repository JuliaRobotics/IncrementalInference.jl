# Showcasing the available analysis tools for the Bayes (Junction) tree.
using IncrementalInference
using DistributedFactorGraphs # For `isSolvable` function.
using RoME # For loading canonical graphs (e.g., Kaess' example).
using Combinatorics # For creating the variable ordering `permutations`.
using Gadfly

latex_fonts = Theme(major_label_font="CMU Serif", major_label_font_size=16pt,
                    minor_label_font="CMU Serif", minor_label_font_size=14pt,
                    key_title_font="CMU Serif", key_title_font_size=12pt,
                    key_label_font="CMU Serif", key_label_font_size=10pt)
Gadfly.push_theme(latex_fonts)

# Get tree for each variable ordering in a factor graph.
fg = loadCanonicalFG_Kaess()
all_trees = getAllTrees(fg)

# scores stores: (tree key ID, nnz, cost fxn 1, cost fxn 2).
unsorted_scores = Vector{Tuple{Int, Float64, Float64, Float64}}()
for key in keys(all_trees)
    e = all_trees[key] # (Bayes tree, var order, nnz)
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

bincnt = 20
Gadfly.plot(x=all_nnzs, y=costs_02,
            Geom.hexbin(xbincount=bincnt, ybincount=bincnt),
            Guide.xlabel("Number of non zeros [int]"),
            Guide.ylabel("Tree cost [cfxn2]"))

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
