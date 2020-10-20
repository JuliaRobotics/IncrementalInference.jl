
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






## =============================================================================
## To be deprecated for consolidation, can revisit after all is one function #675
## =============================================================================

# UP fetch + condition

function getMsgUpChannel(cdat::BayesTreeNodeData)
  Base.depwarn("Deprecated - using Edge Messages for consolidation, #675 - 2", :getMsgUpChannel)
  return cdat.upMsgChannel
end
getMsgUpChannel(cliq::TreeClique) = getMsgUpChannel(getCliqueData(cliq))

function putCliqueMsgUp!(cdat::BayesTreeNodeData, upmsg::LikelihoodMessage)
  Base.depwarn("Deprecated - using take! model", :putCliqueMsgUp!)
  # new replace put! interface
  cdc_ = getMsgUpChannel(cdat)
  if isready(cdc_)
    # first clear an existing value
    take!(cdc_)
  end
  put!(cdc_, upmsg)
  # cdat.upMsg = msg
end

## =============================================================================
## To be deprecated fetch
## =============================================================================


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

  Base.depwarn("Deprecated - using take! model", :prepPutCliqueStatusMsgUp!)
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

# DOWN fetch + condition
function getDwnMsgConsolidated(btnd::BayesTreeNodeData)
  Base.depwarn("Deprecated - using Edge Messages for consolidation, #675 - 2", :getDwnMsgConsolidated)
  btnd.dwnMsgChannel
end
getDwnMsgConsolidated(cliq::TreeClique) = getDwnMsgConsolidated(getCliqueData(cliq))

function fetchDwnMsgConsolidated(btnd::BayesTreeNodeData)
  Base.depwarn("Deprecated - using take! model", :fetchDwnMsgConsolidated)
  fetch(getDwnMsgConsolidated(btnd))
end
fetchDwnMsgConsolidated(cliq::TreeClique) = fetchDwnMsgConsolidated(getCliqueData(cliq))



function putDwnMsgConsolidated!(btnd::BayesTreeNodeData, msg::LikelihoodMessage)
  Base.depwarn("Deprecated - using take! model", :putDwnMsgConsolidated!)
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

Consolidated downward messages generator and Channel sender.

Notes
- Post #459
"""
function prepPutCliqueStatusMsgDwn!(csmc::CliqStateMachineContainer,
                                    status::Symbol=getCliqueStatus(csmc.cliq);
                                    dfg::AbstractDFG=csmc.cliqSubFg,
                                    childSolvDims::Dict{Int,Float64} = Dict{Int,Float64}(),
                                    dwnmsg=prepSetCliqueMsgDownConsolidated!(dfg, csmc.cliq, LikelihoodMessage(status=status, childSolvDims=childSolvDims), csmc.logger, status=status )  )
  #
  Base.depwarn("Deprecated - using take! model", :prepPutCliqueStatusMsgDwn!)

  cd = getCliqueData(csmc.cliq)

  setCliqueStatus!(csmc.cliq, status)

  # NOTE consolidate with upMsgChannel #459
  putDwnMsgConsolidated!(cd, dwnmsg)

  notify(getSolveCondition(csmc.cliq))
  # double notification fixes the problem with hex init (likely due to old lock or something, remove with care)
  sleep(0.1)
  notify(getSolveCondition(csmc.cliq))

  infocsm(csmc, "prepPutCliqueStatusMsgDwn! -- notified status=$(dwnmsg.status), msgs $(collect(keys(dwnmsg.belief))), childSolvDims=$childSolvDims")

  status
end


"""
    $(SIGNATURES)

Return the last up message stored in This `cliq` of the Bayes (Junction) tree.
"""
fetchMsgUpThis(cdat::BayesTreeNodeData) = fetch(getMsgUpChannel(cdat))  # cdat.upMsg    # TODO rename to fetchMsgUp
fetchMsgUpThis(cliql::TreeClique) = fetchMsgUpThis(getCliqueData(cliql))
fetchMsgUpThis(btl::AbstractBayesTree, frontal::Symbol) = fetchMsgUpThis(getClique(btl, frontal))



function fetchMsgsUpChildrenDict( treel::AbstractBayesTree,
                                  cliq::TreeClique )
  #
  msgs = Dict{Int, LikelihoodMessage}()
  for chld in getChildren(treel, cliq)
    msgs[chld.index] = fetchMsgUpThis(chld)
  end

  return msgs
end
fetchMsgsUpChildrenDict( csmc::CliqStateMachineContainer ) = fetchMsgsUpChildrenDict( csmc.tree, csmc.cliq )


"""
    $SIGNATURES

Get and return upward belief messages as stored in child cliques from `treel::AbstractBayesTree`.

Notes
- Use last parameter to select the return format.
- Pull model #674

DevNotes
- Consolidate fetchChildrenStatusUp, getMsgsUpInitChildren
- FIXME update refactor to fetch or take, #855

Related

fetchMsgsUpChildrenDict
"""
function fetchMsgsUpChildren( treel::AbstractBayesTree,
                              cliq::TreeClique,
                              ::Type{TreeBelief} )
  #
  chld = getChildren(treel, cliq)
  retmsgs = Vector{LikelihoodMessage}(undef, length(chld))
  for i in 1:length(chld)
    retmsgs[i] = fetchMsgUpThis(chld[i])
  end
  return retmsgs
end


function fetchMsgsUpChildren( csmc::CliqStateMachineContainer,
                              ::Type{TreeBelief}=TreeBelief )
  #
  # TODO, replace with single channel stored in csmcs or cliques
  fetchMsgsUpChildren(csmc.tree, csmc.cliq, TreeBelief)
end



## TODO Consolidate/Deprecate below


"""
    $SIGNATURES

Fetch (block) caller until child cliques of `cliq::TreeClique` have valid csm status.

Notes:
- Return `::Dict{Symbol}` indicating whether next action that should be taken for each child clique.
- See status options at `getCliqueStatus(..)`.
- Can be called multiple times
"""
function fetchChildrenStatusUp( tree::AbstractBayesTree,
                                cliq::TreeClique,
                                logger=ConsoleLogger() )
  #
  ret = Dict{Int, Symbol}()
  chlr = getChildren(tree, cliq)
  for ch in chlr
      # # FIXME, why are there two steps getting cliq status????
      # chst = getCliqueStatus(ch)  # TODO, remove this
    with_logger(logger) do
      @info "cliq $(cliq.index), child $(ch.index) isready(initUpCh)=$(isready(getMsgUpChannel(ch)))."
    end
    flush(logger.stream)
    # either wait to fetch new result, or report or result
    ret[ch.index] = (fetch(getMsgUpChannel(ch))).status
    with_logger(logger) do
      @info "ret[$(ch.index)]=$(ret[ch.index])."
    end
  end

  return ret
end

# FIXME TEMPORARY CONSOLIDATION FUNCTIONS
# this method adds children and own up msg info to the returning Dict.
# own information is added to capture information from cousins during down init.
function getMsgsUpInitChildren( treel::AbstractBayesTree,
                                cliq::TreeClique,
                                ::Type{TreeBelief};
                                skip::Vector{Int}=Int[])
  #
  chld = getChildren(treel, cliq)
  retmsgs = Dict{Int, LikelihoodMessage}()
  # add possible information that may have come via grandparents from elsewhere in the tree
  if !(cliq.index in skip)
    thismsg = fetchMsgUpThis(cliq)
    retmsgs[cliq.index] = thismsg
  end

  # now add information from each of the child cliques (no longer all stored in prnt i.e. old push #674)
  for ch in chld
    chmsg = fetchMsgUpThis(ch)
    if !(ch.index in skip)
      retmsgs[ch.index] = chmsg
    end
  end
  return retmsgs
end

function getMsgsUpInitChildren( csmc::CliqStateMachineContainer,
                                ::Type{TreeBelief}=TreeBelief;
                                skip::Vector{Int}=Int[] )
  #
  # TODO, replace with single channel stored in csmcs or cliques
  getMsgsUpInitChildren(csmc.tree, csmc.cliq, TreeBelief, skip=skip)
end


