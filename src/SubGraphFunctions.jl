
"""
    $(SIGNATURES)
Specialized subgraph function for cliques to build a deep subgraph copy from the DFG given a list of frontals and separators.
Dev notes:
- TODO Since a clique should already have a list of frontals, seperators, and potentials (factors), this function should just be a light wrapper around copyGraph or buildSubgraph
- TODO Send in clique and then extract frontals, separators and factors
"""
function buildCliqSubgraph!(cliqSubFg::AbstractDFG,
                            dfg::AbstractDFG,
                            frontals::Vector{Symbol},
                            separators::Vector{Symbol};
                            solvable::Int = 0,
                            verbose::Bool=false )


  allvars = union(frontals, separators)
  lenbefore = length(allvars)
  # filter variables by solvable
  solvable != 0 && filter!(fid -> (getSolvable(dfg, fid) >= solvable),  allvars)

  # Potential problem... what are variables doing in the clique if they are not solvable?
  solvable != 0 && lenbefore != length(allvars) && @info("Not all variables are included in subgraph due to solvable $solvable")

  #get list of factors to possibly add, ie. frontal neighbors
  #todo replace with the factor list (potentials) from the clique
  addfac = Symbol[]
  for sym in frontals
    union!(addfac, getNeighbors(dfg,sym))
  end

  allfacs = Symbol[]
  for sym in addfac
    vos = getVariableOrder(dfg, sym)
    if vos ⊆ allvars #only add if not orphaned
      union!(allfacs, [sym])
    end
  end

  # filter factors by solvable
  solvable != 0 && filter!(fid -> (getSolvable(dfg, fid) >= solvable),  allfacs)

  # add all the factors and variables to the new subgraph
  DFG.deepcopyGraph!(cliqSubFg, dfg, allvars, allfacs, verbose=verbose)

  return cliqSubFg
end

function buildCliqSubgraph!(cliqSubFg::AbstractDFG,
                            dfg::AbstractDFG,
                            cliq::TreeClique;
                            solvable::Int = 0,
                            verbose::Bool=false )

  vars = getCliqVarIdsAll(cliq)
  facs = getCliqFactorIdsAll(cliq)
  # Potential problem... what are variables/factors doing in the clique if they are not solvable?
  solvable != 0 && filter!(fid -> (getSolvable(dfg, fid) >= solvable),  vars)
  solvable != 0 && filter!(fid -> (getSolvable(dfg, fid) >= solvable),  facs)

  # fix for issue #681
  # @show ls(cliqSubFg), vars
  if length(intersect(ls(cliqSubFg), vars)) != length(vars)
    DFG.deepcopyGraph!(cliqSubFg, dfg, vars, facs, verbose=verbose)
  end

  return cliqSubFg
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

DevNotes
- TODO review, are all updates atomic?? Then perhaps in-memory only can be reduced to references back to csmc.dfg.
"""
function buildCliqSubgraph(dfg::AbstractDFG,
                           cliq::TreeClique,
                           subfg::InMemoryDFGTypes=InMemDFGType(solverParams=getSolverParams(dfg));
                           solvable::Int=1,
                           verbose::Bool=false )

  #TODO why was solvable hardcoded to 1?
  buildCliqSubgraph!(subfg, dfg, cliq, solvable=solvable, verbose=verbose)
  return subfg
end

function buildCliqSubgraph(fgl::AbstractDFG,
                           treel::AbstractBayesTree,
                           cliqsym::Symbol,
                           subfg::InMemoryDFGTypes=InMemDFGType(solverParams=getSolverParams(fgl));
                           solvable::Int=1,
                           verbose::Bool=false )
  #
  buildCliqSubgraph!(subfg, fgl, getClique(treel, cliqsym), solvable=solvable, verbose=verbose)
  return subfg
end

#


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
                                 updatePPE::Bool=true,
                                 solveKey::Symbol=:default)
  #
  with_logger(logger) do
    @info "transferUpdateSubGraph! -- syms=$syms"
  end

  # transfer specific fields into dest from src
  for var in (x->getVariable(src, x)).(syms)
    # copy not required since a broadcast is used internally
    updateVariableSolverData!(dest, var, solveKey, false, [:val; :bw; :inferdim; :solvedCount]; warn_if_absent=false)
    updatePPE && DFG.updatePPE!(dest, var, solveKey; warn_if_absent=false)
  end

  nothing
end
#     # TODO add with DFG v0.4
#     for sym in syms
#       vari = DFG.getVariable(src, sym)
#       rc = size(getSolverData(vari).val)
#       # TODO -- reduce to DFG functions only
#       pp = getKDE(vari)
#       rc2 = size(getPoints(pp))
#       @info "sym=$sym, mem size of val=$rc and $(rc2)"
#       updateFullVertData!(dest, vari, updatePPE=updatePPE)
#     end
