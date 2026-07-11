
# likely to be deleted at some point

## =============================================================================
## Clique status accessors
## =============================================================================

"""
    $SIGNATURES

Return `::Symbol` status a particular clique is in, with specific regard to solution
or numerical initialization status:
- :needdownmsg
- UPSOLVED
- DOWNSOLVED
- INITIALIZED
- MARGINALIZED
- NULL

Notes:
- `NULL` represents the first uninitialized state of a cliq.
"""
getCliqueStatus(cliqdata::BayesTreeNodeData) = cliqdata.status
getCliqueStatus(cliq::TreeClique) = getCliqueStatus(getCliqueData(cliq))

"""
    $SIGNATURES

Set up initialization or solve status of this `cliq`.
"""
function setCliqueStatus!(cdat::BayesTreeNodeData, status::CliqStatus)
  return cdat.status = status
end
function setCliqueStatus!(cliq::TreeClique, status::CliqStatus)
  return setCliqueStatus!(getCliqueData(cliq), status)
end

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

## ----------------------------------------------------------------------------- 
## DOWN
## ----------------------------------------------------------------------------- 
"""
$SIGNATURES
Get the message channel
"""
getMsgDwnChannel(tree::MetaBayesTree, edge) = MetaGraphs.get_prop(tree.bt, edge, :downMsg)

"""
    $SIGNATURES

Put a belief message on the down tree message channel edge. Blocks until a take! is performed by a different task.
"""
function putBeliefMessageDown!(tree::AbstractBayesTree, edge, beliefMsg::LikelihoodMessage)
  # Blocks until data is available.
  put!(getMsgDwnChannel(tree, edge), beliefMsg)
  return beliefMsg
end

"""
    $SIGNATURES

Remove and return a belief message from the down tree message channel edge. Blocks until data is available.
"""
function takeBeliefMessageDown!(tree::AbstractBayesTree, edge)
  # Blocks until data is available.
  beliefMsg = take!(getMsgDwnChannel(tree, edge))
  return beliefMsg
end

##==============================================================================
## Clique Message Buffers
##==============================================================================
"""
    $SIGNATURES
Get the message buffer that is used to store messages in the clique between state machine steps.
Tx messages are currently used for debugging only and messages are recalculated on each pass. 
"""
function getMessageBuffer(btnd::BayesTreeNodeData)
  return btnd.messages
end
getMessageBuffer(clique::TreeClique) = getCliqueData(clique).messages

"""
    $SIGNATURES
The up message received on the clique.
"""
getMessageUpRx(clique::TreeClique) = getMessageBuffer(clique).upRx
"""
    $SIGNATURES
The down message received on the clique.
"""
getMessageDownRx(clique::TreeClique) = getMessageBuffer(clique).downRx
"""
    $SIGNATURES
Used for debugging only, the up message that was sent by the clique
"""
getMessageUpTx(clique::TreeClique) = getMessageBuffer(clique).upTx
"""
    $SIGNATURES
Used for debugging only, the down message that was sent by the clique
"""
getMessageDownTx(clique::TreeClique) = getMessageBuffer(clique).downTx
