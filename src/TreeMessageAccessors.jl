
export
  getCliqueStatus,
  setCliqueStatus!,
  getSolveCondition

# Reguler accessors
export
  getMsgUpThis,
  fetchDwnMsgsThis,
  fetchMsgDwnInit,
  fetchMsgUpInit

# TODO consolidate
export
  getMsgUpThisInit,
  getMsgUpInitChannel_,
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


getMsgUpChannel(tree::BayesTree, edge) = tree.messages[edge.index].upMsg
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


"""
    $(SIGNATURES)

Set the downward passing message for Bayes (Junction) tree clique `cliql`.
"""
function putMsgDwnThis!(cliql::TreeClique, msgs::LikelihoodMessage)
  getCliqueData(cliql).dwnMsg = msgs
end

function putMsgDwnThis!(csmc::CliqStateMachineContainer, msgs::LikelihoodMessage)
  putMsgDwnThis!(csmc.cliq, msgs)  # NOTE, old, csmc.msgsDown = msgs
end


"""
    $(SIGNATURES)

Return the last down message stored in `cliq` of Bayes (Junction) tree.
"""
fetchMsgDwnThis(cliql::TreeClique) = getCliqueData(cliql).dwnMsg
fetchMsgDwnThis(csmc::CliqStateMachineContainer) = getMsgsDwnThis(csmc.cliq)
fetchMsgDwnThis(btl::AbstractBayesTree, sym::Symbol) = getMsgsDwnThis(getClique(btl, sym))



## =============================================================================
## Consolidating INIT up and down Message Registers/Channels, getters and setters
## =============================================================================


"""
    $(SIGNATURES)

Return the last up message stored in This `cliq` of the Bayes (Junction) tree.
"""
getMsgUpThis(cdat::BayesTreeNodeData) = cdat.upMsg
getMsgUpThis(cliql::TreeClique) = getMsgUpThis(getCliqueData(cliql))
getMsgUpThis(btl::AbstractBayesTree, frontal::Symbol) = getMsgUpThis(getClique(btl, frontal))


"""
    $SIGNATURES

Based on a push model from child cliques that should have already completed their computation.

Dev Notes
- FIXME: Old style -- design has been changed to a Pull model #674
"""
getMsgUpThisInit(cdat::BayesTreeNodeData) = cdat.upMsg # cdat.upInitMsgs
getMsgUpThisInit(cliq::TreeClique) = getMsgUpThisInit(getCliqueData(cliq))

function setCliqueMsgUp!(cdat::BayesTreeNodeData, msg::LikelihoodMessage)
  cdat.upMsg = msg
end

getMsgUpInitChannel_(cdat::BayesTreeNodeData) = cdat.initUpChannel

getMsgDwnThisInit(cdat::BayesTreeNodeData) = cdat.downInitMsg
getMsgDwnInitChannel_(cdat::BayesTreeNodeData) = cdat.initDownChannel

getMsgDwnThisInit(cliq::TreeClique) = getMsgDwnThisInit(getCliqueData(cliq))

getMsgDwnInitChannel_(cliq::TreeClique) = getMsgDwnInitChannel_(getCliqueData(cliq))
fetchMsgDwnInit(cliq::TreeClique) = fetch(getMsgDwnInitChannel_(cliq))

getMsgUpInitChannel_(cliq::TreeClique) = getMsgUpInitChannel_(getCliqueData(cliq))
fetchMsgUpInit(cliq::TreeClique) = fetch(getMsgUpInitChannel_(cliq))


"""
    $(SIGNATURES)

Set the upward passing message for This `cliql` in Bayes (Junction) tree.

Dev Notes
- TODO setUpMsg! should also set inferred dimension
"""
function putMsgUpThis!(cliql::TreeClique,
                       msgs::LikelihoodMessage )
  #
  cd = getCliqueData(cliql)

  # FIXME older interface, likely to be removed at end of #459 and only use upMsgChannel
  setCliqueMsgUp!(cd, msgs)

  # new replace put! interface
  if isready(cd.upMsgChannel)
    # first clear an existing value
    take!(cd.upMsgChannel)
  end
  # insert the new value
  put!(cd.upMsgChannel, msgs)
  nothing
end


function blockMsgDwnUntilStatus(cliq::TreeClique, status::CliqStatus)
  while fetchMsgDwnInit(cliq).status != status
    wait(getSolveCondition(cliq))
  end
  nothing
end

"""
    $SIGNATURES

Blocking call until `cliq` upInit processes has arrived at a result.
"""
function getCliqInitUpResultFromChannel(cliq::TreeClique)
  status = take!(getMsgUpInitChannel_(cliq)).status
  @info "$(current_task()) Clique $(cliq.index), dumping up init status $status"
  return status
end




function putMsgDwnInitStatus!(cliq::TreeClique, status::CliqStatus, logger=ConsoleLogger())
  cdat = getCliqueData(cliq)
  cdc = getMsgDwnInitChannel_(cdat)
    if isready(cdc)
      content = take!(cdc)
      with_logger(logger) do
        @info "dumping stale cliq=$(cliq.index) status message $(content), replacing with $(status)"
      end
    end
  put!(cdc, LikelihoodMessage(status=status))
  notify(getSolveCondition(cliq))
    # FIXME hack to mitigate old race condition
    sleep(0.1)
    notify(getSolveCondition(cliq))

  nothing
end



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
                                   dfg::AbstractDFG=csmc.cliqSubFg)
  #
  # TODO replace with msg channels only

  # construct init's up msg from initialized separator variables
  upinitmsg = prepCliqInitMsgsUp(dfg, csmc.cliq, status)

  # put the init upinitmsg
  putMsgUpThis!(csmc.cliq, upinitmsg )
  # putMsgUpInit!(csmc.cliq, upinitmsg, csmc.logger)

  cond = getSolveCondition(csmc.cliq)
  if getCliqueStatus(csmc.cliq) != status
    infocsm(csmc, "prepPutCliqueStatusMsgUp! -- notify status=$status")
      cdat = getCliqueData(csmc.cliq)
      cdc = getMsgUpInitChannel_(cdat)
        if isready(cdc)
          content = take!(cdc)
        end
      # TODO, lock should not be required in all cases.
      # FIXME reduce hard print and stream flush calls here
      lockUpStatus!(csmc.cliq, csmc.cliq.index, true, csmc.logger, true, "putMsgUpInitStatus!")
      setCliqueStatus!(cdat, status)
      put!(cdc, LikelihoodMessage(status=status))
      notify(cond)
        # FIXME hack to avoid a race condition  -- remove with atomic lock logic upgrade
        sleep(0.1)
      #
      unlockUpStatus!(cdat)
  end
  notify(cond)

  # print a little late
  infocsm(csmc, "8g, doCliqUpsSolveInit. -- postupinitmsg with $(collect(keys(upinitmsg.belief)))")

  # return new up messages in case the user wants to see
  return upinitmsg
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
    retmsgs[i] = getMsgsUpThis(chld[i])
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
function getMsgsUpInitChildren(treel::AbstractBayesTree,
                               cliq::TreeClique,
                               ::Type{TreeBelief},
                               skip::Vector{Int}=Int[])
  #
  chld = getChildren(treel, cliq)
  retmsgs = Dict{Int, LikelihoodMessage}()
  # add possible information that may have come via grandparents from elsewhere in the tree
  thismsg = getMsgUpThisInit(cliq) # TODO change to getmsgUpThis
  retmsgs[cliq.index] = thismsg

  # now add information from each of the child cliques (no longer all stored in prnt i.e. old push #674)
  for ch in chld
    chmsg = getMsgUpThisInit(ch) # TODO change to getmsgUpThis
    if !(ch.index in skip)
      retmsgs[ch.index] = chmsg
    end
  end
  return retmsgs
end

function getMsgsUpInitChildren(csmc::CliqStateMachineContainer,
                               ::Type{TreeBelief}=TreeBelief,
                               skip::Vector{Int}=Int[] )
  #
  # TODO, replace with single channel stored in csmcs or cliques
  getMsgsUpInitChildren(csmc.tree, csmc.cliq, TreeBelief, skip)
end



"""
    $SIGNATURES

Get the latest down message from the parent node (without calculating anything).

Notes
- Different from down initialization messages that do calculate new values -- see `prepCliqInitMsgsDown!`.
- Basically converts function `getDwnMsgs` from `Dict{Symbol,BallTreeDensity}` to `Dict{Symbol,Vector{BallTreeDensity}}`.
"""
function getMsgDwnParent(treel::AbstractBayesTree, cliq::TreeClique)
  downmsgs = IntermediateMultiSiblingMessages()
  for prnt in getParent(treel, cliq)
    for (key, bel) in getDwnMsgs(prnt)
      if !haskey(downmsgs, key)
        downmsgs[key] = IntermediateSiblingMessages()
      end
      # TODO insert true inferred dim
      push!(downmsgs[key], bel)
    end
  end
  return downmsgs
end



#
