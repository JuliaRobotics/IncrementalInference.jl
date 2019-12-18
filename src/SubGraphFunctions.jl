
# TODO -- convert to use add_vertex! instead, since edges type must be made also
function addVerticesSubgraph(fgl::FactorGraph,
    fgseg::FactorGraph,
    vertdict::Dict{Int,Graphs.ExVertex})

    for vert in vertdict
      fgseg.g.vertices[vert[1]] = vert[2]
      if haskey(fgl.v,vert[1])
        fgseg.g.vertices[vert[1]] = vert[2]
        fgseg.IDs[Symbol(vert[2].label)] = vert[1]

        # add edges going in opposite direction
        elr = Graphs.out_edges(vert[2], fgl.g)
        len = length(elr)
        keeprm = trues(len)
        j = 0
        for i in 1:len
          if !haskey(vertdict, elr[i].target.index) # a function node in set, so keep ref
            keeprm[i] = false
            j+=1
          end
        end
        if j < len
          elridx = elr[1].source.index
          fgseg.g.inclist[elridx] = elr[keeprm]
        end
      elseif haskey(fgl.f, vert[1])
        fgseg.f[vert[1]] = vert[2] # adding element to subgraph
        fgseg.fIDs[Symbol(vert[2].label)] = vert[1]
        # get edges associated with function nodes and push edges onto incidence list
        el = Graphs.out_edges(vert[2], fgl.g)
        elidx = el[1].source.index
        fgseg.g.inclist[elidx] = el # okay because treating function nodes only
        fgseg.g.nedges += length(el)
      else
        error("Unknown type factor graph vertex type, something is wrong")
      end
    end
    nothing
end

"""
    $SIGNATURES

Construct a new factor graph object as a subgraph of `fgl::FactorGraph` based on the
variable labels `syms::Vector{Symbols}`.

Notes
- Slighly messy internals, but gets the job done -- some room for performance improvement.

Related

getVariableIds
"""
function buildSubgraphFromLabels(dfg::G,
                                 syms::Vector{Symbol},
                                 destType::Type{<:AbstractDFG}=InMemDFGType;
                                 solvable::Int=0  ) where G <: AbstractDFG

  # data structure for cliq sub graph
  if G <: InMemoryDFGTypes
    #Same type
    cliqSubFg = initfg(G, params=getSolverParams(dfg))
  else
    #Default
    cliqSubFg = initfg(destType, params=getSolverParams(dfg))
  end

  # add a little too many variables (since we need the factors)
  for sym in syms
    DFG.getSubgraphAroundNode(dfg, DFG.getVariable(dfg, sym), 2, false, cliqSubFg, solvable=solvable)
  end

  # remove excessive variables that were copied by neighbors distance 2
  currVars = DFG.getVariableIds(cliqSubFg)
  toDelVars = setdiff(currVars, syms)
  for dv in toDelVars
    # delete any neighboring factors first
    for fc in DFG.lsf(cliqSubFg, dv)
      DFG.deleteFactor!(cliqSubFg, fc)
    end

    # and the variable itself
    DFG.deleteVariable!(cliqSubFg, dv)
  end

  return cliqSubFg
end

#TODO solvable
function buildSubgraphFromLabels!(dfg::AbstractDFG,
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
    if !exists(cliqSubFg,fac) && vos ⊆ allvars   #TODO don't add duplicates to start with
      DFG.addFactor!(cliqSubFg, fac._variableOrderSymbols, deepcopy(fac))
    end
  end

  # remove orphans
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

Transfer contents of `src::FactorGraph` variables `syms::Vector{Symbol}` to `dest::FactorGraph`.

Notes
- Reads, `dest` := `src`, for all `syms`
"""
function transferUpdateSubGraph!(dest::G1,
                                 src::G2,
                                 syms::Vector{Symbol}=union(ls(src)...),
                                 logger=ConsoleLogger()  ) where {G1 <: AbstractDFG, G2 <: AbstractDFG}
  #
  with_logger(logger) do
    @info "transferUpdateSubGraph! -- syms=$syms"

    # TODO add with DFG v0.4
    # DFG.updateGraphSolverData!(src, dest, syms)
    for sym in syms
      vari = DFG.getVariable(src, sym)
      rc = size(solverData(vari).val)
      pp = getKDE(vari)
      rc2 = size(getPoints(pp))
      @info "sym=$sym, mem size of val=$rc and $(rc2)"
      updateFullVertData!(dest, vari, updateMAPest=true)
    end
  end
  nothing
end
