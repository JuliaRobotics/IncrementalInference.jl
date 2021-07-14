"""
    $SIGNATURES

Return true or false depending on whether child cliques are all up solved.
"""
function areCliqChildrenAllUpSolved(treel::AbstractBayesTree,
                                    prnt::TreeClique)::Bool
  #
  for ch in getChildren(treel, prnt)
    if !isCliqUpSolved(ch)
      return false
    end
  end
  return true
end

"""
    $SIGNATURES

Return `true` if any of the children cliques have status `:needdownmsg`.
"""
function doAnyChildrenNeedDwnMsg(children::Vector{TreeClique})::Bool
  for ch in children
    if getCliqueStatus(ch) == :needdownmsg
      return true
    end
  end
  return false
end

function doAnyChildrenNeedDwnMsg(tree::AbstractBayesTree, cliq::TreeClique)::Bool
  doAnyChildrenNeedDwnMsg( getChildren(tree, cliq) )
end

"""
    $SIGNATURES

Return true if has parent with status `:needdownmsg`.
"""
function isCliqParentNeedDownMsg(tree::AbstractBayesTree, cliq::TreeClique, logger=ConsoleLogger())
  prnt = getParent(tree, cliq)
  if length(prnt) == 0
    return false
  end
  prstat = getCliqueStatus(prnt[1])
  with_logger(logger) do
    @info "$(current_task()) Clique $(cliq.index), isCliqParentNeedDownMsg -- parent status: $(prstat)"
  end
  return prstat == :needdownmsg
end

"""
    $SIGNATURES

Wait here if all siblings and the parent status are `:needdownmsg`.
Return true when parent is `INITIALIZED` after all were `:needdownmsg`

Notes
- used for regulating long need down message chains.
- exit strategy is parent becomes status `INITIALIZED`.
"""
function blockCliqSiblingsParentNeedDown( tree::AbstractBayesTree,
                                          cliq::TreeClique,
                                          prnt_::TreeClique; 
                                          logger=ConsoleLogger())
  #
  
  allneeddwn = true
  prstat = getCliqueStatus(prnt_)
  if prstat == :needdownmsg
    for ch in getChildren(tree, prnt_)
      chst = getCliqueStatus(ch)
      if chst != :needdownmsg
        allneeddwn = false
        break;
      end
    end

    if allneeddwn
      # do actual fetch
      prtmsg = fetchDwnMsgConsolidated(prnt_).status
      if prtmsg == INITIALIZED
        return true
      end
    end
  end
  return false
end

"""
    $SIGNATURES

Return a vector of clique index `::Int` for descending order solvable dimensions of each clique.

Notes
- Uses the number of inferable/solvable dimension information in descending order.
- Uses function fetchCliqSolvableDims/getCliqVariableMoreInitDims to retrieved cached values of new solvable/inferable dimensions.

Related

fetchCliqSolvableDims, getCliqVariableMoreInitDims, getSubFgPriorityInitOrder
"""
function getCliqSiblingsPriorityInitOrder(tree::AbstractBayesTree,
                                          prnt::TreeClique,
                                          dwninitmsgs::LikelihoodMessage,
                                          logger=ConsoleLogger() )::Vector{Int}
  #
  sibs = getChildren(tree, prnt)
  len = length(sibs)
  tdims = Vector{Int}(undef, len)
  sidx = Vector{Int}(undef, len)
  for idx in 1:len
    cliqd = getCliqueData(sibs[idx])
    with_logger(logger) do
      @info "$(now()), getCliqSiblingsPriorityInitOrder, sidx=$sidx of $len, $(cliqd.frontalIDs[1]), dwninitmsgs.childSolvDims=$(dwninitmsgs.childSolvDims)"
    end
    flush(logger.stream)
    sidx[idx] = sibs[idx].id
    
        # NEW accumulate the solvableDims for each sibling (#910)
        # FIXME rather determine this using static tree structure and values from dwnmsg (#910)
        @show sidx, dwninitmsgs.childSolvDims
        # @show haskey(dwninitmsgs.childSolvDims,sidx) ? dwninitmsgs.childSolvDims[sidx] : 0

    # FIXME comment out as part of #910
    sidims = fetchCliqSolvableDims(sibs[idx])
    accmSolDims = collect(values(sidims))
    @show tdims[idx] = sum(accmSolDims)
    with_logger(logger) do
      @info "$(now()), prnt=getCliqSiblingsPriorityInitOrder, sidx=$sidx of $len, tdims[idx]=$(tdims[idx])"
    end
    flush(logger.stream)
  end
  p = sortperm(tdims, rev=true)
  with_logger(logger) do
    @info "getCliqSiblingsPriorityInitOrder, done p=$p"
  end
  return sidx[p]
end
