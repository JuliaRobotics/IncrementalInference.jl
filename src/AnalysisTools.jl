"""
    $SIGNATURES

Prunes factor graph to keep up to `upto` number of variables.

Warning: uses functions that are outside of IncrementalInference (e.g.,
`isSolvable()`), so will probably need to place this elsewhere.
"""
function shrinkFactorGraph(fg; upto::Int=6)
  fgs = deepcopy(fg)

  delVars = filter(x->isSolvable(getVariable(fgs, x))==0, ls(fgs))
  todel = setdiff(lsf(fgs, solvable=0), lsf(fgs, solvable=1))
  delFcts = intersect(lsf(fgs), todel)
  allMags = filter(x->:MAGNETOMETER in getTags(getFactor(fgs, x)), lsfPriors(fgs) )
  union!(delFcts, filter(x->length(ls(fgs, x))==0, allMags) )

  union!(delVars, (ls(fgs, r"x\d") |> sortDFG)[upto:end])
  union!(delFcts, map(x->ls(fgs, x), delVars)...)

  map(x->deleteFactor!(fgs, x), delFcts)
  map(x->deleteVariable!(fgs, x), delVars)

  fgs
end


"""
    $SIGNATURES

Deterministically get all trees associated with all possible variable orderings
of `dfg` factor graph. Returns a dictionary with (tree, ordering) tuples.

Warning: factorial number of possibilities, so use carefully!
"""
function getAllTrees(fg::AbstractDFG)
    dfg = loadCanonicalFG_Kaess()
    variables = ls(dfg)
    orderings = permutations(variables) |> collect

    # Dimensionality check to make sure we do not break your computer.
    max_dimension = 11 # something reasonable (11! ~ 40M).
    if length(variables) > 11
        throw(ArgumentError("You crazy! dfg is too big. Factorial explosion."))
    end

    # Produce a tree for each ordering, and store in dictionary.
    all_trees = Dict{Int, Tuple{BayesTree, Vector{Symbol}}}()
    for i in 1:length(orderings)
        all_trees[i] = (resetBuildTreeFromOrder!(fg, orderings[i]), orderings[i])
    end

    return all_trees
end


"""
    $SIGNATURES

Get number of non-zero entries for clique's frontal components. Num of non-zero
matrix entries is just the fully dense upper triangular part of square matrix.
"""
function nnzFrontals(dimension)
    if dimension == 1
        return 1
    else
        # Solved recurrence for n + (n-1) + ... + 2 + 1.
        return (dimension * (dimension + 1)) / 2
    end
end



"""
    $SIGNATURES

Get total number of non-zero entries for a clique. Num of non-zero matrix
entries is the fully dense upper triangular part (frontals) plus the
(frontal x separator)-size rectangle.
"""
function nnzClique(clique)
    frontal_dim = length(getCliqFrontalVarIds(clique))
    separator_dim = length(getCliqSeparatorVarIds(clique))
    return nnzFrontals(frontal_dim) + (frontal_dim * separator_dim)
end


"""
    $SIGNATURES

Get total number of non-zero entries for a Bayes tree. Num of non-zero matrix
entries is the sum of all non-zero entries for each individual clique.
"""
function nnzTree(tree::BayesTree)
  nnzTot = 0
  for (cliqid, cliq) in tree.cliques
    nnzTot += nnzClique(cliq)
  end
  return nnzTot
end
