
export
  getCliqueStatus,
  setCliqueStatus!,
  getSolveCondition

# Reguler accessors
export
  getMsgUpThis,
  fetchDwnMsgsThis,
  fetchMsgDwnInit

# TODO consolidate
export
  getMsgDwnThisInit,
  getMsgDwnInitChannel_

export
  putMsgUpThis!

export
  getMsgDownParent,
  getMsgsUpChildren,
  stackCliqUpMsgsByVariable,
  getCliqDownMsgsAfterDownSolve

# likely to be deleted at some point
export getMsgsUpInitChildren

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
## Clique condition locks
## =============================================================================

# lockUpStatus!(cdat::BayesTreeNodeData, idx::Int=1, verbose::Bool=false, logger=SimpleLogger(stdout), flushLogger::Bool=false, msg::AbstractString="") = put!(cdat.lockUpStatus, idx)
# lockUpStatus!(cliq::TreeClique, idx::Int=1, verbose::Bool=false, logger=SimpleLogger(stdout), flushLogger::Bool=false, msg::AbstractString="") = lockUpStatus!(getCliqueData(cliq), idx, verbose, logger, flushLogger, msg)
function lockUpStatus!(cdat::BayesTreeNodeData,
                       owner::Int=1,
                       idx::Int=1,
                       verbose::Bool=false,
                       logger=SimpleLogger(stdout),
                       flushLogger::Bool=false,
                       msg::AbstractString="")
  #
  with_logger(logger) do
    tt = now()
    verbose ? @info("$(tt) -- Lock $owner, UP from $idx | $msg") : nothing
    verbose && isready(cdat.lockUpStatus) ? @info("$(tt) -- Lock $owner, UP from $idx blocked by $(fetch(cdat.lockUpStatus))") : nothing
  end
  flushLogger ? (flush(logger.stream); sleep(0.01)) : nothing
  put!(cdat.lockUpStatus, idx)
end
lockUpStatus!(cliq::TreeClique, idx::Int=cliq.index, verbose::Bool=false, logger=SimpleLogger(stdout), flushLogger::Bool=false, msg::AbstractString="") = lockUpStatus!(getCliqueData(cliq), cliq.index, idx, verbose, logger, flushLogger, msg)
unlockUpStatus!(cdat::BayesTreeNodeData) = take!(cdat.lockUpStatus)
unlockUpStatus!(cliq::TreeClique) = unlockUpStatus!(getCliqueData(cliq))

function lockDwnStatus!(cdat::BayesTreeNodeData, idx::Int=1; logger=ConsoleLogger())
  with_logger(logger) do
    @info "lockDwnStatus! isready=$(isready(cdat.lockDwnStatus)) with data $(cdat.lockDwnStatus.data)"
  end
  flush(logger.stream)
  stf =  put!(cdat.lockDwnStatus, idx)
  with_logger(logger) do
    @info "lockDwnStatus! DONE: isready=$(isready(cdat.lockDwnStatus)) with data $(cdat.lockDwnStatus.data)"
  end
  flush(logger.stream)
  stf
end
unlockDwnStatus!(cdat::BayesTreeNodeData) = take!(cdat.lockDwnStatus)



## =============================================================================
## Regular up and down Message Registers/Channels, getters and setters
## =============================================================================

## =============================================================================
## Consolidating INIT up and down Message Registers/Channels, getters and setters
## =============================================================================


getMsgUpChannel(tree::BayesTree, edge) = tree.messages[edge.index].upMsg  # FIXME convert to Channel only
getMsgDwnChannel(tree::BayesTree, edge) = tree.messages[edge.index].downMsg

getMsgUpChannel(tree::MetaBayesTree, edge) = MetaGraphs.get_prop(tree.bt, edge, :upMsg)
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

Remove and return belief message from the up tree message channel edge. Blocks until data is available.
"""
function takeBeliefMessageUp!(tree::AbstractBayesTree, edge)
  # Blocks until data is available.
  beliefMsg = take!(getMsgUpChannel(tree, edge))
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

Put a belief message on the up tree message channel `edge`. Blocks until a take! is performed by a different task.
"""
function putBeliefMessageUp!(tree::BayesTree, edge, beliefMsg::LikelihoodMessage)
  # Blocks until data is available.
  put!(getMsgUpChannel(tree, edge), beliefMsg)
  return beliefMsg
end


getMsgUpChannel(cdat::BayesTreeNodeData) = cdat.upMsgChannel
getMsgUpChannel(cliq::TreeClique) = getMsgUpChannel(getCliqueData(cliq))



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
    $(SIGNATURES)

Return the last up message stored in This `cliq` of the Bayes (Junction) tree.
"""
getMsgUpThis(cdat::BayesTreeNodeData) = fetch(getMsgUpChannel(cdat))  # cdat.upMsg    # TODO rename to fetchMsgUp
getMsgUpThis(cliql::TreeClique) = getMsgUpThis(getCliqueData(cliql))
getMsgUpThis(btl::AbstractBayesTree, frontal::Symbol) = getMsgUpThis(getClique(btl, frontal))



"""
    $SIGNATURES

Notify of new up status and message.

DevNotes
- Major part of #459 consolidation effort.
- FIXME require consolidation
  - Perhaps deprecate putMsgUpInitStatus! as separate function?
"""
function prepPutCliqueStatusMsgUp!(csmc::CliqStateMachineContainer,
                                   status::Symbol=getCliqueStatus(csmc.cliq);
                                   dfg::AbstractDFG=csmc.cliqSubFg,
                                   upmsg=prepCliqInitMsgsUp(dfg, csmc.cliq, status)  )
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

  infocsm(csmc, "prepPutCliqueStatusMsgUp! -- notified status=$status with msg keys $(collect(keys(upmsg.belief)))")

  # return new up messages in case the user wants to see
  return upmsg
end


## =============================================================================
## Family message getters and setters
## =============================================================================


"""
    $SIGNATURES

Get and return upward belief messages as stored in child cliques from `treel::AbstractBayesTree`.

Notes
- Use last parameter to select the return format.
- Pull model #674

DevNotes
- Consolidate two versions getMsgsUpChildren
"""
function getMsgsUpChildren(treel::AbstractBayesTree,
                           cliq::TreeClique,
                           ::Type{TreeBelief} )
  #
  chld = getChildren(treel, cliq)
  retmsgs = Vector{LikelihoodMessage}(undef, length(chld))
  for i in 1:length(chld)
    retmsgs[i] = getMsgUpThis(chld[i])
  end
  return retmsgs
end


function getMsgsUpChildren(csmc::CliqStateMachineContainer,
                           ::Type{TreeBelief}=TreeBelief )
  #
  # TODO, replace with single channel stored in csmcs or cliques
  getMsgsUpChildren(csmc.tree, csmc.cliq, TreeBelief)
end

# FIXME TEMPORARY CONSOLIDATION FUNCTIONS
# this method adds children and own up msg info to the returning Dict.
# own information is added to capture information from cousins during down init.
function getMsgsUpInitChildren(treel::AbstractBayesTree,
                               cliq::TreeClique,
                               ::Type{TreeBelief};
                               skip::Vector{Int}=Int[])
  #
  chld = getChildren(treel, cliq)
  retmsgs = Dict{Int, LikelihoodMessage}()
  # add possible information that may have come via grandparents from elsewhere in the tree
  if !(cliq.index in skip)
    thismsg = getMsgUpThis(cliq)
    retmsgs[cliq.index] = thismsg
  end

  # now add information from each of the child cliques (no longer all stored in prnt i.e. old push #674)
  for ch in chld
    chmsg = getMsgUpThis(ch)
    if !(ch.index in skip)
      retmsgs[ch.index] = chmsg
    end
  end
  return retmsgs
end

function getMsgsUpInitChildren(csmc::CliqStateMachineContainer,
                               ::Type{TreeBelief}=TreeBelief;
                               skip::Vector{Int}=Int[] )
  #
  # TODO, replace with single channel stored in csmcs or cliques
  getMsgsUpInitChildren(csmc.tree, csmc.cliq, TreeBelief, skip=skip)
end




#
