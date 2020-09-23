
# these functions and their use in the code need to be consolidated to conclude #459

# TODO consolidate
export
  fetchDwnMsgsThis,
  getMsgDwnThisInit,
  getMsgDownParent


## ============================================================================
## Consolidated dwn message container
## ============================================================================

getDwnMsgConsolidated(btnd::BayesTreeNodeData) = btnd.dwnMsgChannel
getDwnMsgConsolidated(cliq::TreeClique) = getDwnMsgConsolidated(getCliqueData(cliq))

function fetchDwnMsgConsolidated(btnd::BayesTreeNodeData)
  fetch(getDwnMsgConsolidated(btnd))
end
fetchDwnMsgConsolidated(cliq::TreeClique) = fetchDwnMsgConsolidated(getCliqueData(cliq))

function putDwnMsgConsolidated!(btnd::BayesTreeNodeData, msg::LikelihoodMessage)
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
                                    dwnmsg=getSetDownMessagesComplete!(dfg, csmc.cliq, LikelihoodMessage(), csmc.logger, status=status )  )
  #
  cd = getCliqueData(csmc.cliq)

  setCliqueStatus!(csmc.cliq, status)

  # NOTE consolidate with upMsgChannel #459
  putDwnMsgConsolidated!(cd, dwnmsg)

  notify(getSolveCondition(csmc.cliq))
  # took ~40 hours to figure out that a double norification fixes the problem with hex init
  sleep(0.1)
  notify(getSolveCondition(csmc.cliq))

  infocsm(csmc, "prepPutCliqueStatusMsgDwn! -- notified status=$(dwnmsg.status) with msg keys $(collect(keys(dwnmsg.belief)))")

  status
end



## ============================================================================
## more channels, MUST BE REMOVED
## ============================================================================


# FIXME DEPRECATE TO RATHER USE getDwnMsgConsolidated
getMsgDwnChannel(tree::BayesTree, edge) = tree.messages[edge.index].downMsg
getMsgDwnChannel(tree::MetaBayesTree, edge) = MetaGraphs.get_prop(tree.bt, edge, :downMsg)

"""
    $SIGNATURES

Remove and return a belief message from the down tree message channel edge. Blocks until data is available.
"""
function takeBeliefMessageDown!(tree::BayesTree, edge)
  # Blocks until data is available.
  beliefMsg = take!(getMsgDwnChannel(tree, edge))
  return beliefMsg
end


"""
    $SIGNATURES

Put a belief message on the down tree message channel edge. Blocks until a take! is performed by a different task.
"""
function putBeliefMessageDown!(tree::BayesTree, edge, beliefMsg::LikelihoodMessage)
  # Blocks until data is available.
  put!(getMsgDwnChannel(tree, edge), beliefMsg)
  return beliefMsg
end





"""
    $SIGNATURES

Calculate new and then set the down messages for a clique in Bayes (Junction) tree.
"""
function getSetDownMessagesComplete!( subfg::AbstractDFG,
                                      cliq::TreeClique,
                                      prntDwnMsgs::LikelihoodMessage,
                                      logger=ConsoleLogger();
                                      status::CliqStatus=getCliqueStatus(cliq)  )
  #
  allvars = getCliqVarIdsAll(cliq)
  allprntkeys = collect(keys(prntDwnMsgs.belief))
  passkeys = intersect(allvars, setdiff(allprntkeys,ls(subfg)))
  remainkeys = setdiff(allvars, passkeys)
  newDwnMsgs = LikelihoodMessage(status=status)

  # some msgs are just pass through from parent
  for pk in passkeys
    newDwnMsgs.belief[pk] = prntDwnMsgs.belief[pk]
  end

  # other messages must be extracted from subfg
  for mk in remainkeys
    setVari = getVariable(subfg, mk)
    if isInitialized(setVari)
      newDwnMsgs.belief[mk] = TreeBelief(setVari)
    end
  end

  # set the downward keys
  with_logger(logger) do
    @info "cliq $(cliq.index), getSetDownMessagesComplete!, allkeys=$(allvars), passkeys=$(passkeys), msgkeys=$(collect(keys(newDwnMsgs.belief)))"
  end

  return newDwnMsgs
end




## =============================================================================
## Atomic messaging during init -- TODO deprecated
## =============================================================================




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





## =============================================================
