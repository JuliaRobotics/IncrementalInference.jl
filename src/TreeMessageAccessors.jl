
export
  getCliqStatus,
  setCliqStatus!,
  getSolveCondition

# Reguler accessors
export
  fetchMsgUpThis,
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
  putMsgUpThis!,
  putMsgUpInit!,
  putMsgUpInitStatus!

export
  getMsgDownParent,
  getMsgsUpChildren,
  stackCliqUpMsgsByVariable,
  getCliqDownMsgsAfterDownSolve

# likely to be deleted at some point
export getMsgsUpChildrenInitDict

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
getCliqStatus(cliqdata::BayesTreeNodeData) = cliqdata.initialized
getCliqStatus(cliq::TreeClique) = getCliqStatus(getCliqueData(cliq))
getCliqStatusUp(cliq::TreeClique) = getCliqStatus(cliq)

"""
    $SIGNATURES

Set up initialization or solve status of this `cliq`.
"""
function setCliqStatus!(cliq::TreeClique, status::Symbol)
  getCliqueData(cliq).initialized = status
end



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

Set the upward passing message for This `cliql` in Bayes (Junction) tree.

Dev Notes
- TODO setUpMsg! should also set inferred dimension
"""
function putMsgUpThis!(cliql::TreeClique, msgs::LikelihoodMessage)
  # TODO change into a replace put!
  cd = getCliqueData(cliql)

  # older interface
  cd.upMsg = msgs

  # new interface
  if isready(cd.upMsgChannel)
    # first clear an existing value
    take!(cd.upMsgChannel)
  end
  # insert the new value
  put!(cd.upMsgChannel, msgs)
  nothing
end


"""
    $(SIGNATURES)

Return the last up message stored in This `cliq` of the Bayes (Junction) tree.
"""
fetchMsgUpThis(cliql::TreeClique) = getCliqueData(cliql).upMsg
fetchMsgUpThis(btl::AbstractBayesTree, frontal::Symbol) = getUpMsgs(getClique(btl, frontal))


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
    $SIGNATURES

Based on a push model from child cliques that should have already completed their computation.

Dev Notes
- FIXME: Old style -- design has been changed to a Pull model #674
"""
getMsgUpThisInit(cdat::BayesTreeNodeData) = cdat.upInitMsgs
getMsgUpInitChannel_(cdat::BayesTreeNodeData) = cdat.initUpChannel

getMsgDwnThisInit(cdat::BayesTreeNodeData) = cdat.downInitMsg
getMsgDwnInitChannel_(cdat::BayesTreeNodeData) = cdat.initDownChannel

getMsgUpThisInit(cliq::TreeClique) = getMsgUpThisInit(getCliqueData(cliq))
getMsgDwnThisInit(cliq::TreeClique) = getMsgDwnThisInit(getCliqueData(cliq))

getMsgDwnInitChannel_(cliq::TreeClique) = getMsgDwnInitChannel_(getCliqueData(cliq))
fetchMsgDwnInit(cliq::TreeClique) = fetch(getMsgDwnInitChannel_(cliq))

getMsgUpInitChannel_(cliq::TreeClique) = getMsgUpInitChannel_(getCliqueData(cliq))
fetchMsgUpInit(cliq::TreeClique) = fetch(getMsgUpInitChannel_(cliq))


function setMsgUpThisInitDict!(cdat::BayesTreeNodeData, idx, msg::LikelihoodMessage)
  getMsgUpThisInit(cdat) = msg
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


"""
    $SIGNATURES

Set cliques up init msgs.

DevNotes
- ORIGINALLY PART OF PUSH MODEL #674, MUST BE UPDATED TO PULL.
  -- Likely problem for siblings wanting to have notified parent
    -- Notifications might have to remain on parent while msgs are stored in each' own clique
- TODO, must be consolidated
"""
function putMsgUpInit!(cliq::TreeClique,
                       childid::Int,
                       msg::LikelihoodMessage,
                       logger=SimpleLogger(stdout))
  #
  cd = getCliqueData(cliq)
  soco = getSolveCondition(cliq)
  # FIXME, locks should not be required in all cases
  lockUpStatus!(cliq, cliq.index, true, logger, true, "putMsgUpInit!") # TODO XX
  # FIXME, consolidation required, convert to Pull model #674
  setMsgUpThisInitDict!(cd, childid, msg)
  # TODO simplify and fix need for repeat
  # notify cliq condition that there was a change
  notify(soco)
  #hack for mitigating deadlocks, in case a user was not already waiting, but waiting on lock instead
  sleep(0.1)
  notify(soco)
  unlockUpStatus!(cd)
  nothing
end


function putMsgUpInitStatus!(cliq::TreeClique, status::CliqStatus, logger=SimpleLogger(stdout))
  cdat = getCliqueData(cliq)
  cdc = getMsgUpInitChannel_(cdat)
  cond = getSolveCondition(cliq)
    if isready(cdc)
      content = take!(cdc)
    end
  # FIXME, lock should not be required in all cases.
  lockUpStatus!(cliq, cliq.index, true, logger, true, "putMsgUpInitStatus!") # TODO XX
  cdat.initialized = status
  put!(cdc, LikelihoodMessage(status=status))
  notify(cond)
    # FIXME hack to avoid a race condition  -- remove with atomic lock logic upgrade
    sleep(0.1)
    notify(cond) # getSolveCondition(cliq)
  #
  unlockUpStatus!(cdat)
  nothing
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
function getMsgsUpChildrenInitDict(treel::AbstractBayesTree,
                                   cliq::TreeClique,
                                   ::Type{TreeBelief},
                                   skip::Vector{Int}=Int[])
  #
  chld = getChildren(treel, cliq)
  retmsgs = Dict{Int, LikelihoodMessage}()
  # add possible information that may have come via grandparents from elsewhere in the tree
  thismsg = getMsgUpThisInit(cliq)
  retmsgs[cliq.index] = thismsg
  # @assert length(thismsg) <= 1 "getMsgUpThisInit must contain this clique local info only."
  # for (ke, va) in thismsg
    # retmsgs[ke] = va
  # end

  # now add information from each of the child cliques (no longer all stored in prnt i.e. old push #674)
  # retmsgs = Vector{LikelihoodMessage}(undef, length(chld))
  for ch in chld
    # @show cliq.index, ch.index, skip, collect(keys(getMsgUpThisInit(ch)))
    chmsg = getMsgUpThisInit(ch)
    # @assert !(length(chmsg) == 1 && !haskey(chmsg, ch.index)) "getMsgUpThisInit must contain only local clique messages."
    # if haskey(chmsg, ch.index) # FIXME, this should not be required, since it wasnt before
    if !(ch.index in skip)
    # if length(chmsg) == 1 && !(ch.index in skip)
      retmsgs[ch.index] = chmsg # [ch.index] # getMsgUpThisInit(ch) # TODO X
    end
    # end
  end
  return retmsgs
end
function getMsgsUpChildrenInitDict(csmc::CliqStateMachineContainer,
                                   ::Type{TreeBelief}=TreeBelief,
                                   skip::Vector{Int}=Int[] )
  #
  # TODO, replace with single channel stored in csmcs or cliques
  getMsgsUpChildrenInitDict(csmc.tree, csmc.cliq, TreeBelief, skip)
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
