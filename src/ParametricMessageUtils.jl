

"""
    $SIGNATURES

Remove and return a belief message from the down tree message channel edge. Blocks until data is available.
"""
function takeBeliefMessageDown!(tree::BayesTree, edge)
  # Blocks until data is available.
  beliefMsg = take!(tree.messages[edge.index].downMsg)
  return beliefMsg
end

function takeBeliefMessageDown!(tree::MetaBayesTree, edge)
  # Blocks until data is available.
  beliefMsg = take!(MetaGraphs.get_prop(tree.bt, edge, :downMsg))
  return beliefMsg
end


"""
    $SIGNATURES

Remove and return belief message from the up tree message channel edge. Blocks until data is available.
"""
function takeBeliefMessageUp!(tree::BayesTree, edge)
  # Blocks until data is available.
  beliefMsg = take!(tree.messages[edge.index].upMsg)
  return beliefMsg
end

function takeBeliefMessageUp!(tree::MetaBayesTree, edge)
  # Blocks until data is available.
  beliefMsg = take!(MetaGraphs.get_prop(tree.bt, edge, :upMsg))
  return beliefMsg
end


"""
    $SIGNATURES

Put a belief message on the down tree message channel edge. Blocks until a take! is performed by a different task.
"""
function putBeliefMessageDown!(tree::BayesTree, edge, beliefMsg::LikelihoodMessage)
  # Blocks until data is available.
  put!(tree.messages[edge.index].downMsg, beliefMsg)
  return beliefMsg
end

function putBeliefMessageDown!(tree::MetaBayesTree, edge, beliefMsg::LikelihoodMessage)
  # Blocks until data is available.
  put!(MetaGraphs.get_prop(tree.bt, edge, :downMsg), beliefMsg)
  return beliefMsg
end


"""
    $SIGNATURES

Put a belief message on the up tree message channel `edge`. Blocks until a take! is performed by a different task.
"""
function putBeliefMessageUp!(tree::BayesTree, edge, beliefMsg::LikelihoodMessage)
  # Blocks until data is available.
  put!(tree.messages[edge.index].upMsg, beliefMsg)
  return beliefMsg
end

function putBeliefMessageUp!(tree::MetaBayesTree, edge, beliefMsg::LikelihoodMessage)
  # Blocks until data is available.
  put!(MetaGraphs.get_prop(tree.bt, edge, :upMsg), beliefMsg)
  return beliefMsg
end
