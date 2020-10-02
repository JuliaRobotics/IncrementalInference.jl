
# these functions need to be consolidated to single get and set API


## ============================================================================
## Consolidated dwn message GETS -- where msgs are stored
## ============================================================================

getDwnMsgConsolidated(btnd::BayesTreeNodeData) = btnd.dwnMsgChannel
getDwnMsgConsolidated(cliq::TreeClique) = getDwnMsgConsolidated(getCliqueData(cliq))

# FIXME CONSOLIDATE TO SINGLE CHANNEL STORAGE LOCATION
getDwnMsgConsolidated(tree::BayesTree, edge) = tree.messages[edge.index].downMsg
getDwnMsgConsolidated(tree::MetaBayesTree, edge) = MetaGraphs.get_prop(tree.bt, edge, :downMsg)


## ============================================================================
## Dwn message FETCH / TAKE  -- TODO consolidate to single function
## ============================================================================


function fetchDwnMsgConsolidated(btnd::BayesTreeNodeData)
  fetch(getDwnMsgConsolidated(btnd))
end
fetchDwnMsgConsolidated(cliq::TreeClique) = fetchDwnMsgConsolidated(getCliqueData(cliq))


"""
    $SIGNATURES

Remove and return a belief message from the down tree message channel edge. Blocks until data is available.
"""
function takeBeliefMessageDown!(tree::BayesTree, edge)
  # Blocks until data is available.
  beliefMsg = take!(getDwnMsgConsolidated(tree, edge))
  return beliefMsg
end



## =============================================================================
## Down message SETS / PUTS -- TODO consolidate to a single function (#855)
## =============================================================================


"""
    $SIGNATURES

Put a belief message on the down tree message channel edge. Blocks until a take! is performed by a different task.
"""
function putBeliefMessageDown!(tree::BayesTree, edge, beliefMsg::LikelihoodMessage)
  # Blocks until data is available.
  put!(getDwnMsgConsolidated(tree, edge), beliefMsg)
  return beliefMsg
end


function putDwnMsgConsolidated!(btnd::BayesTreeNodeData, msg::LikelihoodMessage)
  # need to get the current solvableDims
  dmc = getDwnMsgConsolidated(btnd)
  if isready(dmc)
    take!(dmc)
  end
  put!(dmc, msg)
end
putDwnMsgConsolidated!(cliq::TreeClique, msg::LikelihoodMessage) = putDwnMsgConsolidated!(getCliqueData(cliq), msg)


"""
    $SIGNATURES

THIS IS ONE OF THE FAVORITES FOR POST CONSOLIDATED DOWNWARD MESSAGES.
"""
function prepPutCliqueStatusMsgDwn!(csmc::CliqStateMachineContainer,
                                    status::Symbol=getCliqueStatus(csmc.cliq);
                                    dfg::AbstractDFG=csmc.cliqSubFg,
                                    dwnmsg=prepSetCliqueMsgDownConsolidated!(dfg, csmc.cliq, LikelihoodMessage(status=status), csmc.logger, status=status )  )
  #
  cd = getCliqueData(csmc.cliq)

  setCliqueStatus!(csmc.cliq, status)

  # NOTE consolidate with upMsgChannel #459
  putDwnMsgConsolidated!(cd, dwnmsg)

  notify(getSolveCondition(csmc.cliq))
  # double notification fixes the problem with hex init (likely due to old lock or something, remove with care)
  sleep(0.1)
  notify(getSolveCondition(csmc.cliq))

  infocsm(csmc, "prepPutCliqueStatusMsgDwn! -- notified status=$(dwnmsg.status) with msg keys $(collect(keys(dwnmsg.belief)))")

  status
end




# FIXME, consolidate into single down msg event API
function notifyCliqDownInitStatus!( cliq::TreeClique,
                                    status::Symbol;
                                    logger=ConsoleLogger() )
  #
  cdat = getCliqueData(cliq)
    with_logger(logger) do
    @info "$(now()) $(current_task()), cliq=$(cliq.index), notifyCliqDownInitStatus! -- pre-lock, new $(cdat.initialized)-->$(status)"
  end

  # take lock for atomic transaction
  # lockDwnStatus!(cdat, cliq.index, logger=logger)

  setCliqueStatus!(cdat, status)

  # TODO, should this not send the beliefs aswell??
  msg = LikelihoodMessage(status=status)
  putDwnMsgConsolidated!(cliq, msg)
  # putMsgDwnInitStatus!(cliq, status, logger, msg)
  

  # unlock for others to proceed
  # unlockDwnStatus!(cdat)
  with_logger(logger) do
    @info "$(now()), cliq=$(cliq.index), notifyCliqDownInitStatus! -- unlocked, $(getCliqueStatus(cliq))"
  end

  # flush(logger.stream)

  nothing
end


# FIXME figure out if this can be consolidated and deprecated
function dwnPrepOutMsg( fg::AbstractDFG,
                        cliq::TreeClique,
                        dwnMsgs::Array{LikelihoodMessage,1},
                        d::Dict{Symbol, T},
                        logger=ConsoleLogger()) where T
  # pack all downcoming conditionals in a dictionary too.
  with_logger(logger) do
    if cliq.index != 1 #TODO there may be more than one root
      @info "Dwn msg keys $(keys(dwnMsgs[1].belief))"
      @info "fg vars $(ls(fg))"
    end # ignore root, now incoming dwn msg
  end
  m = LikelihoodMessage()
  i = 0
  for vid in getCliqueData(cliq).frontalIDs
    m.belief[vid] = deepcopy(d[vid]) # TODO -- not sure if deepcopy is required
  end
  for cvid in getCliqueData(cliq).separatorIDs
    i+=1
    # TODO -- convert to points only since kde replace by rkhs in future
    m.belief[cvid] = deepcopy(dwnMsgs[1].belief[cvid]) # TODO -- maybe this can just be a union(,)
  end
  return m
end




## =============================================================
