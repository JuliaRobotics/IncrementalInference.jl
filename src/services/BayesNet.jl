
"""
    $SIGNATURES

Determine the variable ordering used to construct both the Bayes Net and Bayes/Junction/Elimination tree.

Notes
- Heuristic method -- equivalent to QR or Cholesky.
- Are using Blas `QR` function to extract variable ordering.
- **NOT USING SUITE SPARSE** -- which would requires commercial license.
- For now `A::Array{<:Number,2}` as a dense matrix.
- Columns of `A` are system variables, rows are factors (without differentiating between partial or full factor).
- default is to use `solvable=1` and ignore factors and variables that might be used for dead reckoning or similar.

Future
- TODO: `A` should be sparse data structure (when we exceed 10'000 var dims)
- TODO: Incidence matrix is rectagular and adjacency is the square.
"""
function getEliminationOrder(
  dfg::AbstractDFG;
  ordering::Symbol = :qr,
  solvable::Int = 1,
  constraints::Vector{Symbol} = Symbol[],
)
  #
  @assert 0 == length(constraints) || ordering == :ccolamd "Must use ordering=:ccolamd when trying to use constraints"
  # Get the sparse adjacency matrix, variable, and factor labels
  adjMat, permuteds, permutedsf = DFG.getBiadjacencyMatrix(dfg; solvable = solvable)
  # adjMat, permuteds, permutedsf = DFG.getAdjacencyMatrixSparse(dfg, solvable=solvable)

  # Create dense adjacency matrix

  p = Int[]
  if ordering == :chol
    # hack for dense matrix....
    A = adjMat
    p = cholesky(Matrix(A'A), Val(true)).piv
    @warn "check that cholesky ordering is not reversed -- basically how much fill in (separator size) are you seeing???  Long skinny chains in tree is bad."
  elseif ordering == :qr
    # hack for dense matrix....
    A = Array(adjMat)
    # this is the default
    q, r, p = qr(A, (v"1.7" <= VERSION ? ColumnNorm() : Val(true)))
    p .= p |> reverse
  elseif ordering == :ccolamd
    cons = zeros(Int, length(adjMat.colptr) - 1)
    cons[findall(x -> x in constraints, permuteds)] .= 1
    p = _ccolamd(adjMat, cons)
      # cons = zeros(SuiteSparse_long, length(adjMat.colptr) - 1)
      # cons[findall(x -> x in constraints, permuteds)] .= 1
      # p = Ccolamd.ccolamd(adjMat, cons)
    @warn "Integration via AMD.ccolamd under development and replaces pre-Julia 1.9 direct ccall approach." maxlog=5
  elseif ordering == :mcs
    # maximum cardinality search
    A = adjMat
    p, _ = CliqueTrees.permutation(A'A; alg=CliqueTrees.MCS())
  elseif ordering == :rcm
    # reverse Cuthill-Mckee
    A = adjMat
    p, _ = CliqueTrees.permutation(A'A; alg=CliqueTrees.RCM())
  elseif ordering == :mmd
    # multiple minimum degree
    A = adjMat
    p, _ = CliqueTrees.permutation(A'A; alg=CliqueTrees.MMD())
  else
    @error("getEliminationOrder -- cannot do the requested ordering $(ordering)")
  end

  # Return the variable ordering that we should use for the Bayes map
  # reverse order checked in #475 and #499
  return permuteds[p]
end

# lets create all the vertices first and then deal with the elimination variables thereafter
function addBayesNetVerts!(dfg::AbstractDFG, elimOrder::Array{Symbol, 1})
  #
  for pId in elimOrder
    vert = DFG.getVariable(dfg, pId)
    if  getVariableState(vert).BayesNetVertID === nothing ||
        getVariableState(vert).BayesNetVertID == :_null # Special serialization case of nothing
      @debug "[AddBayesNetVerts] Assigning $pId.data.BayesNetVertID = $pId"
      getVariableState(vert).BayesNetVertID = pId
    else
      @warn "addBayesNetVerts -- Something is wrong, variable '$pId' should not have an existing Bayes net reference to '$(getVariableState(vert).BayesNetVertID)'"
    end
  end
end

function addConditional!(dfg::AbstractDFG, vertId::Symbol, Si::Vector{Symbol})
  #
  bnv = DFG.getVariable(dfg, vertId)
  bnvd = getVariableState(bnv)
  bnvd.separator = Si
  for s in Si
    push!(bnvd.BayesNetOutVertIDs, s)
  end
  return nothing
end

function addChainRuleMarginal!(dfg::AbstractDFG, Si::Vector{Symbol})
  #

  lbls = String[]
  genmarg = GenericMarginal()
  Xi = map(v -> DFG.getVariable(dfg, v), Si)
  # @info "adding marginal to"
  # for x in Xi
  #   @info "x.index=",x.index
  # end
  addFactor!(dfg, Xi, genmarg; graphinit = false, suppressChecks = true)
  return nothing
end

function rmVarFromMarg(dfg::AbstractDFG, fromvert::DFGVariable, gm::Vector{DFGFactor})
  #

  @debug " - Removing $(fromvert.label)"
  for m in gm
    @debug "Looking at $(m.label)"
    for n in listNeighbors(dfg, m) #x1, x2
      if n == getLabel(fromvert) # n.label ==? x1
        @debug "   - Breaking link $(m.label)->$(fromvert.label)..."
        @debug "     - Original links: $(DFG.ls(dfg, m))"
        remvars = setdiff(DFG.ls(dfg, m), [fromvert.label])
        @debug "     - New links: $remvars"

        DFG.deleteFactor!(dfg, m) # Remove it
        if length(remvars) > 0
          @debug "$(m.label) still has links to other variables, readding it back..."
          addFactor!(
            dfg,
            remvars,
            _getCCW(m).usrfnc!;
            graphinit = false,
            suppressChecks = true,
          )
        else
          @debug "$(m.label) doesn't have any other links, not adding it back..."
        end
      end
    end
    # Added back in chain rule.
    if DFG.exists(dfg, m) && length(listNeighbors(dfg, m)) <= 1
      @warn "removing vertex id=$(m.label)"
      DFG.deleteFactor!(dfg, m)
    end
  end
  return nothing
end

function buildBayesNet!(dfg::AbstractDFG, elimorder::Vector{Symbol}; solvable::Int = 1)
  #
  # addBayesNetVerts!(dfg, elimorder)
  for v in elimorder
    @debug """ 
                Eliminating $(v)
                ===============
          """
    # which variable are we eliminating

    # all factors adjacent to this variable
    fi = Symbol[]
    Si = Symbol[]
    gm = DFGFactor[]

    vert = DFG.getVariable(dfg, v)
    for fctId in listNeighbors(dfg, vert; solvable = solvable)
      fct = DFG.getFactor(dfg, fctId)
      if (DFG.getFactorState(fct).eliminated != true)
        push!(fi, fctId)
        for sepNode in listNeighbors(dfg, fct; solvable = solvable)
          # TODO -- validate !(sepNode.index in Si) vs. older !(sepNode in Si)
          if sepNode != v && !(sepNode in Si) # Symbol comparison!
            push!(Si, sepNode)
          end
        end
        DFG.getFactorState(fct).eliminated = true
      end

      if typeof(_getCCW(fct)) == CommonConvWrapper{GenericMarginal}
        push!(gm, fct)
      end
    end

    if v != elimorder[end]
      addConditional!(dfg, v, Si)
      # not yet inserting the new prior p(Si) back into the factor graph
    end

    # mark variable
    getVariableState(vert).eliminated = true

    # TODO -- remove links from current vertex to any marginals
    rmVarFromMarg(dfg, vert, gm)

    #add marginal on remaining variables... ? f(xyz) = f(x | yz) f(yz)
    # new function between all Si (round the outside, right the outside)
    length(Si) > 0 && addChainRuleMarginal!(dfg, Si)
  end
  return nothing
end
