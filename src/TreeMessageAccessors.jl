
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
## Message channel put/take! + buffer message accessors
## =============================================================================

## ----------------------------------------------------------------------------- 
## UP
## ----------------------------------------------------------------------------- 
"""
$SIGNATURES
Get the message channel
"""
getMsgUpChannel(tree::BayesTree, edge) = tree.messageChannels[edge.index].upMsg
getMsgUpChannel(tree::MetaBayesTree, edge) = MetaGraphs.get_prop(tree.bt, edge, :upMsg)

"""
$SIGNATURES

Put a belief message on the up tree message channel `edge`. Blocks until a take! is performed by a different task.
"""
function putBeliefMessageUp!(tree::AbstractBayesTree, edge, beliefMsg::LikelihoodMessage)
  # Blocks until data is available.
  put!(getMsgUpChannel(tree, edge), beliefMsg)
  return beliefMsg
end

"""
$SIGNATURES

Remove and return belief message from the up tree message channel edge. Blocks until data is available.
"""
function takeBeliefMessageUp!(tree::AbstractBayesTree, edge)
  # Blocks until data is available.
  beliefMsg = take!(getMsgUpChannel(tree, edge))
  return beliefMsg
end

# @deprecate putBeliefMessageUp!(tree::AbstractBayesTree, edge, beliefMsg::LikelihoodMessage) putMessageUp!(tree, edge, beliefMsg)
# @deprecate takeBeliefMessageUp!(tree::AbstractBayesTree, edge) takeMessageUp!(tree, edge)
# @deprecate putBeliefMessageDown!(tree::BayesTree, edge, beliefMsg::LikelihoodMessage) putMessageDown!(tree, edge, beliefMsg)
# @deprecate takeBeliefMessageDown!(tree::BayesTree, edge) takeMessageDown!(tree, edge)

## ----------------------------------------------------------------------------- 
## DOWN
## ----------------------------------------------------------------------------- 
"""
$SIGNATURES
Get the message channel
"""
getMsgDwnChannel(tree::BayesTree, edge) = tree.messageChannels[edge.index].downMsg
getMsgDwnChannel(tree::MetaBayesTree, edge) = MetaGraphs.get_prop(tree.bt, edge, :downMsg)

@deprecate getDwnMsgConsolidated(tree::AbstractBayesTree, edge) getMsgDwnChannel(tree, edge)

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

Remove and return a belief message from the down tree message channel edge. Blocks until data is available.
"""
function takeBeliefMessageDown!(tree::BayesTree, edge)
  # Blocks until data is available.
  beliefMsg = take!(getMsgDwnChannel(tree, edge))
  return beliefMsg
end


##==============================================================================
## Clique Message Buffers
##==============================================================================
function getMessageBuffer(btnd::BayesTreeNodeData)
  btnd.messages
end
getMessageBuffer(clique::TreeClique) = getCliqueData(clique).messages

# getMessageUpRx(clique::TreeClique) = getMessageBuffer(clique).upRx
# getMessageDownRx(clique::TreeClique) = getMessageBuffer(clique).downRx

# getMessageUpTx(clique::TreeClique) = getMessageBuffer(clique).upTx
# getMessageDownTx(clique::TreeClique) = getMessageBuffer(clique).downTx


