#TODO JT can be removed, used as sanity check
function removeSeparatorPriorsFromSubgraph!(cliqSubFg::AbstractDFG, cliq::TreeClique)
  cliqSeparatorVarIds = getCliqSeparatorVarIds(cliq)
  priorIds = Symbol[]
  for v in cliqSeparatorVarIds
    facs = getNeighbors(cliqSubFg, v)
    for f in facs
      isprior = length(getFactor(cliqSubFg, f)._variableOrderSymbols) == 1
      isprior && push!(priorIds, f)
      isprior && DFG.deleteFactor!(cliqSubFg, f)
    end
  end
  return priorIds
end

# TODO JT solvable
function buildCliqSubgraph!(dfg::AbstractDFG,
                            cliqSubFg::AbstractDFG,
                            frontals::Vector{Symbol},
                            separators::Vector{Symbol};
                            solvable::Int=0)

  for sym in separators
    DFG.addVariable!(cliqSubFg, deepcopy(DFG.getVariable(dfg, sym)))
  end

  addfac = Symbol[]
  for sym in frontals
    DFG.addVariable!(cliqSubFg, deepcopy(DFG.getVariable(dfg, sym)))
    append!(addfac, getNeighbors(dfg,sym))
  end

  allvars = ls(cliqSubFg)
  for sym in addfac
    fac = DFG.getFactor(dfg, sym)
    vos = fac._variableOrderSymbols
    if !exists(cliqSubFg,fac) && vos ⊆ allvars   #duplicates not added to start with
      DFG.addFactor!(cliqSubFg, fac._variableOrderSymbols, deepcopy(fac))
    end
  end

  #TODO remove, just as a sanity check to see if there is any orphans to remove
  for fct in DFG.getFactors(cliqSubFg)
    # delete any neighboring factors first
    if length(getNeighbors(cliqSubFg, fct)) != length(fct._variableOrderSymbols)
      DFG.deleteFactor!(cliqSubFg, fc)
      @error "deleteFactor! this should not happen"
    end
  end

  return cliqSubFg
end



"""
    $SIGNATURES
Transfer contents of `src::AbstractDFG` variables `syms::Vector{Symbol}` to `dest::AbstractDFG`.
Notes
- Reads, `dest` := `src`, for all `syms`
"""
function transferUpdateSubGraph!(dest::AbstractDFG,
                                 src::AbstractDFG,
                                 syms::Vector{Symbol}=union(ls(src)...),
                                 logger=ConsoleLogger();
                                 updatePPE::Bool=true  )
  #
  with_logger(logger) do
    @info "transferUpdateSubGraph! -- syms=$syms"

    # TODO add with DFG v0.4
    # DFG.updateGraphSolverData!(src, dest, syms)
    for sym in syms
      vari = DFG.getVariable(src, sym)
      rc = size(solverData(vari).val)
      # TODO -- reduce to DFG functions only
      pp = getKDE(vari)
      rc2 = size(getPoints(pp))
      @info "sym=$sym, mem size of val=$rc and $(rc2)"
      updateFullVertData!(dest, vari, updatePPE=updatePPE)
    end
  end
  nothing
end



"""
    $SIGNATURES

Build a new subgraph from `fgl<:AbstractDFG` containing all variables and factors
associated with `cliq`.  Additionally add the upward message prior factors as
needed for belief propagation (inference).

Notes
- `cliqsym::Symbol` defines the cliq where variable appears as a frontal variable.
- `varsym::Symbol` defaults to the cliq frontal variable definition but can in case a
  separator variable is required instead.
"""
function buildCliqSubgraphDown(fgl::AbstractDFG, treel::AbstractBayesTree, cliqsym::Symbol, varsym::Symbol=cliqsym)
  @warn "Obsolete, buildCliqSubGraph*() is no longer in use"
  # build a subgraph copy of clique
  cliq = whichCliq(treel, cliqsym)
  syms = getCliqAllVarIds(cliq)
  subfg = buildSubgraphFromLabels(fgl,syms)

  # add upward messages to subgraph
  msgs = getCliqParentMsgDown(treel, cliq)
  addMsgFactors!(subfg, msgs)
  return subfg
end
