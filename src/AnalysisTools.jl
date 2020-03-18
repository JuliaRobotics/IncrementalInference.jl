export shrinkFactorGraph, getAllTrees, nnzFrontals, nnzClique, nnzTree, nnzSqrtInfoMatrix, getTreeCost_01, getTreeCost_02

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
of `dfg` factor graph. Returns a dictionary with (tree, ordering, nnz) tuples.

Warning: factorial number of possibilities, so use carefully!
"""
function getAllTrees(fg::AbstractDFG)
    # dfg = generateCanonicalFG_Kaess(graphinit=false)
    variables = ls(fg)
    orderings = permutations(variables) |> collect

    # Dimensionality check to make sure we do not break your computer.
    max_dimension = 11 # something reasonable (11! ~ 40M).
    if length(variables) > 11
        throw(ArgumentError("You crazy! dfg is too big. Factorial explosion."))
    end

    # Produce a tree for each ordering, and store in dictionary.
    all_trees = Dict{Int, Tuple{BayesTree, Vector{Symbol}, Float64}}()
    for i in 1:length(orderings)
        tree = resetBuildTreeFromOrder!(fg, orderings[i])
        nnz = nnzTree(tree)
        all_trees[i] = (tree, orderings[i], nnz)
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


"""
    $SIGNATURES

Get total number of non-zero entries for a factor graph's upper triangular
square root information matrix, i.e., R matrix in A = Q[R 0]^T, using the QR's
factorization algorithm variable ordering.
"""
function nnzSqrtInfoMatrix(A::Matrix)
  q,r,p = qr(A, Val(true))
  r .= abs.(r)
  nz = 1e-5 .< r
  r[nz] .= 1
  sum(nz)
end


"""
    $SIGNATURES

Simple cost function for ranking the structure of a Bayes tree. Weighting:
    cost = (max tree depth) * (max clique dimension)^alpha
"""
function getTreeCost_01(tree::BayesTree; alpha::Float64=1.0 )
  cliqs = tree.cliques |> values |> collect
  maxdepth = map(x->getCliqDepth(tree, x)+1, cliqs) |> maximum
  maxdim = length.(map(x->getCliqVarIdsAll(x), cliqs)) |> maximum

  return maxdepth * (maxdim^alpha)
end


"""
    $SIGNATURES

Cost function for ranking the structure of a Bayes tree, putting and emphasis on
wider but shallower trees by penalizing the average number of siblings.
Weighting:
    cost = 1/(total num of child / num of parents) *
                           (max tree depth) * (max clique dimension)^alpha
"""
function getTreeCost_02(tree::BayesTree; alpha::Float64=1.0)
  # Frontal and number of children.
  ARR = Tuple{Symbol, Int}[]
  for (cliqid, vex) in tree.cliques
    afrtl = getCliqFrontalVarIds(tree.cliques[cliqid])[1]
    numch = length(getChildren(tree, afrtl))
    push!(ARR, (afrtl, numch))
  end

  numParents = filter(x->0<x[2], ARR) |> length
  totalNumChildren = (x->x[2]).(ARR) |> sum

  return getTreeCost_01(tree, alpha=alpha)/(totalNumChildren/numParents)
end
