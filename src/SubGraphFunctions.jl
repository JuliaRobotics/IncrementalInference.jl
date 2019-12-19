
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
