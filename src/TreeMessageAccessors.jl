
export
  getCliqueStatus,
  setCliqueStatus!,
  getSolveCondition

export
  stackCliqUpMsgsByVariable,
  getCliqDownMsgsAfterDownSolve

# likely to be deleted at some point

## =============================================================================
## Clique status accessors
## =============================================================================

"""
    $SIGNATURES

Return `::Symbol` status a particular clique is in, with specific regard to solution
or numerical initialization status:
- :needdownmsg
- :upsolved
- :downsolved
- :initialized
- :marginalized
- :null

Notes:
- `:null` represents the first uninitialized state of a cliq.
"""
getCliqueStatus(cliqdata::BayesTreeNodeData) = cliqdata.initialized
getCliqueStatus(cliq::TreeClique) = getCliqueStatus(getCliqueData(cliq))

"""
    $SIGNATURES

Set up initialization or solve status of this `cliq`.
"""
function setCliqueStatus!(cdat::BayesTreeNodeData, status::Symbol)
  cdat.initialized = status
end
setCliqueStatus!(cliq::TreeClique, status::Symbol) = setCliqueStatus!(getCliqueData(cliq), status)



## =============================================================================
## Regular up and down Message Registers/Channels, getters and setters
## =============================================================================

## =============================================================================
## Consolidating INIT up and down Message Registers/Channels, getters and setters
## =============================================================================

# FIXME DEPRECATE to use only MsgUpThis register
getMsgUpChannel(tree::BayesTree, edge) = tree.messages[edge.index].upMsg
getMsgUpChannel(tree::MetaBayesTree, edge) = MetaGraphs.get_prop(tree.bt, edge, :upMsg)

# Consolidate to only use this
getMsgUpChannel(cdat::BayesTreeNodeData) = cdat.upMsgChannel
getMsgUpChannel(cliq::TreeClique) = getMsgUpChannel(getCliqueData(cliq))



"""
    $SIGNATURES

Put a belief message on the up tree message channel `edge`. Blocks until a take! is performed by a different task.
"""
function putBeliefMessageUp!(tree::AbstractBayesTree, edge, beliefMsg::LikelihoodMessage)
  # Blocks until data is available.
  put!(getMsgUpChannel(tree, edge), beliefMsg)
  return beliefMsg
end


function putCliqueMsgUp!(cdat::BayesTreeNodeData, upmsg::LikelihoodMessage)
  # new replace put! interface
  cdc_ = getMsgUpChannel(cdat)
  if isready(cdc_)
    # first clear an existing value
    take!(cdc_)
  end
  put!(cdc_, upmsg)
  # cdat.upMsg = msg
end




"""
    $SIGNATURES

Notify of new up status and message.

DevNotes
- Major part of #459 consolidation effort.
- FIXME require consolidation
  - Perhaps deprecate putMsgUpInitStatus! as separate function?
"""
function prepPutCliqueStatusMsgUp!( csmc::CliqStateMachineContainer,
                                    status::Symbol=getCliqueStatus(csmc.cliq);
                                    dfg::AbstractDFG=csmc.cliqSubFg,
                                    upmsg=prepCliqueMsgUpConsolidated(dfg, csmc.cliq, status, logger=csmc.logger)  )
  #
  # TODO replace with msg channels only

  # put the init upmsg
  cd = getCliqueData(csmc.cliq)

  setCliqueStatus!(csmc.cliq, status)

  # NOTE consolidate with upMsgChannel #459
  putCliqueMsgUp!(cd, upmsg)
    # TODO remove as part of putCliqueMsgUp!
    # new replace put! interface
    # cdc_ = getMsgUpChannel(cd)
    # if isready(cdc_)
    #   # first clear an existing value
    #   take!(cdc_)
    # end
    # put!(cdc_, upmsg)

  notify(getSolveCondition(csmc.cliq))
  # took ~40 hours to figure out that a double norification fixes the problem with hex init
  sleep(0.1)
  notify(getSolveCondition(csmc.cliq))

  # also notify parent as part of upinit (PUSH NOTIFY PARENT), received at `slowWhileInit_StateMachine`
  prnt = getParent(csmc.tree, csmc.cliq)
  0 < length(prnt) ? notify(getSolveCondition(prnt[1])) : nothing

  infocsm(csmc, "prepPutCliqueStatusMsgUp! -- notified status=$(upmsg.status) with msg keys $(collect(keys(upmsg.belief)))")

  # return new up messages in case the user wants to see
  return upmsg
end




#
