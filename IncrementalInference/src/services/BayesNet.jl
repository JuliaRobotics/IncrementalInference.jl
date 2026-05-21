
"""
    $SIGNATURES

Determine the variable ordering used to construct both the Bayes Net and Bayes/Junction/Elimination tree.

Notes
- Heuristic method -- equivalent to QR or Cholesky.
- Are using Blas `QR` function to extract variable ordering.
- **NOT USING SUITE SPARSE** -- which would requires commercial license.
- Columns of `A` are system variables, rows are factors (without differentiating between partial or full factor).
- default is to use `solvable=1` and ignore factors and variables that might be used for dead reckoning or similar.
"""
function getEliminationOrder(
  dfg::AbstractDFG;
  ordering::Symbol = :qr,
  solvable::Int = 1,
  constraints::Vector{Symbol} = Symbol[],
)
  # Get the sparse adjacency matrix, variable, and factor labels
  adjMat, permuteds, permutedsf = DFG.getBiadjacencyMatrix(dfg; solvable = solvable)

  # get constraint indices
  clique = findall(∈(constraints), permuteds)

  if ordering == :ccolamd
    # ccolamd handles constraints internally
    cons = zeros(Int, length(adjMat.colptr) - 1)
    cons[clique] .= 1
    p = _ccolamd(adjMat, cons)
    @warn "Integration via AMD.ccolamd under development and replaces pre-Julia 1.9 direct ccall approach." maxlog=5
  else
    S = adjMat' * adjMat

    # force contraint indices to be a clique
    S[clique, clique] .= 1

    if ordering == :chol
      p = cholesky(Matrix(S), Val(true)).piv
      @warn "check that cholesky ordering is not reversed -- basically how much fill in (separator size) are you seeing???  Long skinny chains in tree is bad."
    elseif ordering == :qr
      A = Array(adjMat)
      q, r, p = qr(A, (v"1.7" <= VERSION ? ColumnNorm() : Val(true)))
      reverse!(p)
    elseif ordering == :mcs
      p, _ = CliqueTrees.permutation(S; alg=CliqueTrees.MCS())
    elseif ordering == :rcm
      p, _ = CliqueTrees.permutation(S; alg=CliqueTrees.RCM())
    elseif ordering == :mmd
      p, _ = CliqueTrees.permutation(S; alg=CliqueTrees.MMD())
    else
      @error("getEliminationOrder -- cannot do the requested ordering $(ordering)")
    end

    # move constraints to end of ordering
    if !isempty(clique)
      p, _ = CliqueTrees.permutation(S; alg=CliqueTrees.CompositeRotations(clique, p))
    end
  end

  # Return the variable ordering that we should use for the Bayes map
  # reverse order checked in #475 and #499
  return permuteds[p]
end

function addConditional!(dfg::AbstractDFG, vertId::Symbol, Si::Vector{Symbol})
  #
  bnv = DFG.getVariable(dfg, vertId)
  bnvd = getState(bnv, :default)
  bnvd.separator = Si
  # for s in Si
  #   push!(bnvd.BayesNetOutVertIDs, s)
  # end
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

function rmVarFromMarg(dfg::AbstractDFG, fromvert::VariableCompute, gm::Vector{FactorCompute})
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
    gm = FactorCompute[]

    vert = DFG.getVariable(dfg, v)
    for fctId in listNeighbors(dfg, vert; whereSolvable = >=(solvable))
      fct = DFG.getFactor(dfg, fctId)
      if (fct.state.eliminated != true)
        push!(fi, fctId)
        for sepNode in listNeighbors(dfg, fct; whereSolvable = >=(solvable))
          # TODO -- validate !(sepNode.index in Si) vs. older !(sepNode in Si)
          if sepNode != v && !(sepNode in Si) # Symbol comparison!
            push!(Si, sepNode)
          end
        end
        fct.state.eliminated = true
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
    # getState(vert, :default).eliminated = true #TODO e8d remove? looks unused

    # TODO -- remove links from current vertex to any marginals
    rmVarFromMarg(dfg, vert, gm)

    #add marginal on remaining variables... ? f(xyz) = f(x | yz) f(yz)
    # new function between all Si (round the outside, right the outside)
    length(Si) > 0 && addChainRuleMarginal!(dfg, Si)
  end
  return nothing
end
