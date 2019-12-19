
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

    # TODO add with DFG v0.4
    # DFG.updateGraphSolverData!(src, dest, syms)
    for sym in syms
      vari = DFG.getVariable(src, sym)
      rc = size(solverData(vari).val)
      # TODO -- reduce to DFG functions only
      pp = getKDE(vari)
      rc2 = size(getPoints(pp))
      @info "sym=$sym, mem size of val=$rc and $(rc2)"
      updateFullVertData!(dest, vari, updateMAPest=true)
    end
  end
  nothing
end


# """
#     $SIGNATURES
#
# Construct a new factor graph object as a subgraph of `fgl::FactorGraph` based on the
# variable labels `syms::Vector{Symbols}`.
#
# Notes
# - Slighly messy internals, but gets the job done -- some room for performance improvement.
#
# Related
#
# getVariableIds
# """
# function buildSubgraphFromLabels!(dfg::AbstractDFG,
#                                   syms::Vector{Symbol},
#                                   destType::Type{<:AbstractDFG}=InMemDFGType;
#                                   solvable::Int=0  )
#
#   # data structure for cliq sub graph
#   if G <: InMemoryDFGTypes
#     #Same type
#     cliqSubFg = initfg(G, params=getSolverParams(dfg))
#   else
#     #Default
#     cliqSubFg = initfg(destType, params=getSolverParams(dfg))
#   end
#
#   # add a little too many variables (since we need the factors)
#   for sym in syms
#     DFG.getSubgraphAroundNode(dfg, DFG.getVariable(dfg, sym), 2, false, cliqSubFg, solvable=solvable)
#   end
#
#   # remove excessive variables that were copied by neighbors distance 2
#   currVars = DFG.getVariableIds(cliqSubFg)
#   toDelVars = setdiff(currVars, syms)
#   for dv in toDelVars
#     # delete any neighboring factors first
#     for fc in DFG.lsf(cliqSubFg, dv)
#       DFG.deleteFactor!(cliqSubFg, fc)
#     end
#
#     # and the variable itself
#     DFG.deleteVariable!(cliqSubFg, dv)
#   end
#
#   return cliqSubFg
# end

#TODO solvable
# function buildSubgraphFromLabels!(dfg::AbstractDFG,
#                                   cliqSubFg::AbstractDFG,
#                                   frontals::Vector{Symbol},
#                                   separators::Vector{Symbol};
#                                   solvable::Int=0)
#
#   for sym in separators
#     DFG.addVariable!(cliqSubFg, deepcopy(DFG.getVariable(dfg, sym)))
#   end
#
#   addfac = Symbol[]
#   for sym in frontals
#     DFG.addVariable!(cliqSubFg, deepcopy(DFG.getVariable(dfg, sym)))
#     append!(addfac, getNeighbors(dfg,sym))
#   end
#
#   allvars = ls(cliqSubFg)
#   for sym in addfac
#     fac = DFG.getFactor(dfg, sym)
#     vos = fac._variableOrderSymbols
#     if !exists(cliqSubFg,fac) && vos ⊆ allvars   #TODO don't add duplicates to start with
#       DFG.addFactor!(cliqSubFg, fac._variableOrderSymbols, deepcopy(fac))
#     end
#   end
#
#   # remove orphans
#   for fct in DFG.getFactors(cliqSubFg)
#     # delete any neighboring factors first
#     if length(getNeighbors(cliqSubFg, fct)) != length(fct._variableOrderSymbols)
#       DFG.deleteFactor!(cliqSubFg, fc)
#       @error "deleteFactor! this should not happen"
#     end
#   end
#
#   return cliqSubFg
# end
