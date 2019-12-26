
"""
    $SIGNATURES

Transfer contents of `src::AbstractDFG` variables `syms::Vector{Symbol}` to `dest::AbstractDFG`.

Notes
- Reads, `dest` := `src`, for all `syms`
"""
function transferUpdateSubGraph!(dest::AbstractDFG,
                                 src::AbstractDFG,
                                 syms::Vector{Symbol}=union(ls(src)...),
                                 logger=ConsoleLogger()  )
  #
  with_logger(logger) do
    @info "transferUpdateSubGraph! -- syms=$syms"
  end
  
  # pass workload up to DFG
  DFG.updateGraphSolverData!(src, dest, syms)
  
  nothing
end
# TODO add with DFG v0.4
# for sym in syms
#    vari = DFG.getVariable(src, sym)
#    rc = size(solverData(vari).val)
#    # TODO -- reduce to DFG functions only
#    pp = getKDE(vari)
#    rc2 = size(getPoints(pp))
#    @info "sym=$sym, mem size of val=$rc and $(rc2)"
#    updateFullVertData!(dest, vari, updateMAPest=true)
# end



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
function buildCliqSubgraphDown(fgl::AbstractDFG, treel::BayesTree, cliqsym::Symbol, varsym::Symbol=cliqsym)
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
