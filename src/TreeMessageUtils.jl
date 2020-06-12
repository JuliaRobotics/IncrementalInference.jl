# init utils for tree based inference


## =============================================================================
# helper functions to add tree messages to subgraphs
## =============================================================================


"""
    $SIGNATURES

Modify the `subfg::AbstractDFG` to include `msgs` as priors that are used
during clique inference.

Notes
- May be used initialization or inference, in both upward and downward directions.
- `msgs` are identified by variable label `::Symbol`, and my consist of multiple beliefs.
- Message sets from different cliques are identified by clique id `::Int`.

Related

`deleteMsgFactors!`
"""
function addMsgFactors!(subfg::AbstractDFG,
                        msgs::LikelihoodMessage)::Vector{DFGFactor}
  # add messages as priors to this sub factor graph
  msgfcts = DFGFactor[]
  svars = DFG.listVariables(subfg)
  for (msym, dm) in msgs.belief
    if msym in svars
      # manifold information is contained in the factor graph DFGVariable object
      fc = addFactor!(subfg, [msym],
              MsgPrior(manikde!(dm.val, dm.bw[:,1], getManifolds(dm.softtype)), dm.inferdim), graphinit=false)
      push!(msgfcts, fc)
    end
  end
  return msgfcts
end

function addMsgFactors!(subfg::AbstractDFG,
                        msgs::Dict{Symbol, Vector{Tuple{BallTreeDensity, Float64}}} )
      # msgs::
  # add messages as priors to this sub factor graph
  msgfcts = DFGFactor[]
  svars = DFG.listVariables(subfg)
  for (msym, dms) in msgs
    for dm in dms
      if msym in svars
        # TODO should be on manifold prior, not just generic euclidean prior -- okay since variable on manifold, but not for long term
        fc = addFactor!(subfg, [msym], MsgPrior(dm[1], dm[2]), graphinit=false)
        push!(msgfcts, fc)
      end
    end
  end
  return msgfcts
end

function addMsgFactors!(subfg::AbstractDFG,
                        allmsgs::Dict{Int,LikelihoodMessage} )
  #
  allfcts = DFGFactor[]
  for (cliqid, msgs) in allmsgs
    # do each dict in array separately
    newfcts = addMsgFactors!(subfg, msgs)
    union!( allfcts, newfcts )
  end
  return allfcts
end



# TODO move addMsgFactors_Parametric! here



"""
    $SIGNATURES

Delete from the subgraph`::AbstractDFG` prior belief `msgs` that could/would be used
during clique inference.

Related

`addMsgFactors!`
"""
function deleteMsgFactors!(subfg::AbstractDFG,
                           fcts::Vector{DFGFactor} )
  #
  for fc in fcts
    deleteFactor!(subfg, fc.label)
  end
end



## =============================================================================
# tree belief message accessors for from parameteric development
## =============================================================================



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



## =============================================================================
# tree belief message accessors for nonparameteric in this section
## =============================================================================


"""
    $SIGNATURES

Based on a push model from child cliques that should have already completed their computation.
"""
getCliqInitUpMsgs(cliq::TreeClique)::Dict{Int, LikelihoodMessage} = getCliqueData(cliq).upInitMsgs

getInitDownMsg(cliq::TreeClique)::LikelihoodMessage = getCliqueData(cliq).downInitMsg

"""
    $SIGNATURES

Set cliques up init msgs.
"""
function setCliqUpInitMsgs!(cliq::TreeClique, childid::Int, msg::LikelihoodMessage)
  cd = getCliqueData(cliq)
  soco = getSolveCondition(cliq)
  lockUpStatus!(cd)
  cd.upInitMsgs[childid] = msg
  # notify cliq condition that there was a change
  notify(soco)
  unlockUpStatus!(cd)
  #hack for mitigating deadlocks, in case a user was not already waiting, but waiting on lock instead
  sleep(0.1)
  notify(soco)
  nothing
end


"""
    $(SIGNATURES)

Set the upward passing message for Bayes (Junction) tree clique `cliql`.

Dev Notes
- TODO setUpMsg! should also set inferred dimension
"""
function setUpMsg!(cliql::TreeClique, msgs::LikelihoodMessage)
  # TODO change into a replace put!
  getCliqueData(cliql).upMsg = msgs
end

# NOTE decided not to store messages in CSMC, but closer to Tree instead (likely on edges)
# function setUpMsg!(csmc::CliqStateMachineContainer, cliqid::Int, msgs::LikelihoodMessage)
#   csmc.msgsUp[cliqid] = msgs
# end
# getUpMsgs(csmc::CliqStateMachineContainer) = csmc.msgsUp



"""
    $(SIGNATURES)

Return the last up message stored in `cliq` of Bayes (Junction) tree.
"""
getUpMsgs(cliql::TreeClique) = getCliqueData(cliql).upMsg
getUpMsgs(btl::AbstractBayesTree, sym::Symbol) = getUpMsgs(getCliq(btl, sym))


"""
    $(SIGNATURES)

Set the downward passing message for Bayes (Junction) tree clique `cliql`.
"""
function setDwnMsg!(csmc::CliqStateMachineContainer, msgs::LikelihoodMessage)
  csmc.msgsDown = msgs
end

function setDwnMsg!(cliql::TreeClique, msgs::LikelihoodMessage)
  getCliqueData(cliql).dwnMsg = msgs
end


"""
    $(SIGNATURES)

Return the last down message stored in `cliq` of Bayes (Junction) tree.
"""
getDwnMsgs(csmc::CliqStateMachineContainer) = csmc.msgsDown
getDwnMsgs(cliql::TreeClique) = getCliqueData(cliql).dwnMsg
getDwnMsgs(btl::AbstractBayesTree, sym::Symbol) = getDwnMsgs(getCliq(btl, sym))




"""
    $SIGNATURES

Get and return upward belief messages as stored in child cliques from `treel::AbstractBayesTree`.

Notes
- Use last parameter to select the return format.
- Pull model #674

DevNotes
- Consolidate two versions getCliqChildMsgsUp
"""
function getCliqChildMsgsUp(fg_::AbstractDFG,
                            treel::AbstractBayesTree,
                            cliq::TreeClique,
                            ::Type{TreeBelief} )
  #
  childmsgs = LikelihoodMessage[]
  for child in getChildren(treel, cliq)
    nbpchild = LikelihoodMessage()
    for (key, bel) in getUpMsgs(child).belief
      # manis = getManifolds(fg_, key)
      # inferdim = getVariableInferredDim(fg_, key)
      dcBel = deepcopy(bel)
      nbpchild.belief[key] = TreeBelief(dcBel.val, dcBel.bw, dcBel.inferdim, getSofttype(getVariable(fg_, key)))
    end
    push!(childmsgs, nbpchild)
  end
  return childmsgs
end

function getCliqChildMsgsUp(treel::AbstractBayesTree,
                            cliq::TreeClique,
                            ::Type{BallTreeDensity})
  #
  childmsgs = IntermediateMultiSiblingMessages()
  for child in getChildren(treel, cliq)
    for (key, bel) in getUpMsgs(child).belief
      # id = fg_.IDs[key]
      # manis = getManifolds(fg_, id)
      if !haskey(childmsgs, key)
        childmsgs[key] = IntermediateSiblingMessages()
      end
      push!(childmsgs[key], bel )
    end
  end
  return childmsgs
end

function getCliqChildMsgsUp(csmc::CliqStateMachineContainer,
                            ::Type{TreeBelief} )
  #
  getCliqChildMsgsUp(csmc.cliqSubFg, csmc.tree, csmc.cliq, TreeBelief)
end

function getCliqChildMsgsUp(csmc::CliqStateMachineContainer,
                            ::Type{BallTreeDensity})
  #
  getCliqChildMsgsUp(csmc.tree, csmc.cliq, BallTreeDensity)
end

"""
    $SIGNATURES

Get the latest down message from the parent node (without calculating anything).

Notes
- Different from down initialization messages that do calculate new values -- see `prepCliqInitMsgsDown!`.
- Basically converts function `getDwnMsgs` from `Dict{Symbol,BallTreeDensity}` to `Dict{Symbol,Vector{BallTreeDensity}}`.
"""
function getCliqParentMsgDown(treel::AbstractBayesTree, cliq::TreeClique)
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



"""
    $SIGNATURES

Blocking call until `cliq` upInit processes has arrived at a result.
"""
function getCliqInitUpResultFromChannel(cliq::TreeClique)
  status = take!(getCliqueData(cliq).initUpChannel)
  @info "$(current_task()) Clique $(cliq.index), dumping initUpChannel status, $status"
  return status
end




lockUpStatus!(cdat::BayesTreeNodeData, idx::Int=1) = put!(cdat.lockUpStatus, idx)
lockUpStatus!(cliq::TreeClique, idx::Int=1) = lockUpStatus!(getCliqueData(cliq), idx)
unlockUpStatus!(cdat::BayesTreeNodeData) = take!(cdat.lockUpStatus)
unlockUpStatus!(cliq::TreeClique) = unlockUpStatus!(getCliqueData(cliq))

function lockDwnStatus!(cdat::BayesTreeNodeData, idx::Int=1; logger=ConsoleLogger())
  with_logger(logger) do
    @info "lockDwnStatus! isready=$(isready(cdat.lockDwnStatus)) with data $(cdat.lockDwnStatus.data)"
  end
    flush(logger.stream)
  # if isready(cdat.lockDwnStatus)
  #   with_logger(logger) do
  #     @info "lockDwnStatus! is locked with $(cdat.lockDwnStatus.data)"
  #   end
  #   flush(logger.stream)
  # end
  stf =  put!(cdat.lockDwnStatus, idx)
  with_logger(logger) do
    @info "lockDwnStatus! DONE: isready=$(isready(cdat.lockDwnStatus)) with data $(cdat.lockDwnStatus.data)"
  end
  flush(logger.stream)
  stf
end
unlockDwnStatus!(cdat::BayesTreeNodeData) = take!(cdat.lockDwnStatus)

"""
    $SIGNATURES

Update clique status and notify of the change

Notes
- Assumes users will lock the status state before getting status until after decision whether to update status.
- If so, only unlock after status and condition has been updated.

Dev Notes
- Should be made an atomic transaction
"""
function notifyCliqUpInitStatus!(cliq::TreeClique, status::Symbol; logger=ConsoleLogger())
  cd = getCliqueData(cliq)
  with_logger(logger) do
    tt = split(string(now()), 'T')[end]
    @info "$(tt) $(current_task()), cliq=$(cliq.index), notifyCliqUpInitStatus! -- pre-lock, $(cd.initialized)-->$(status)"
  end
  flush(logger.stream)

  ## TODO only notify if not data structure is not locked by other user (can then remove the hack)
  # Wait until lock can be aquired
  lockUpStatus!(cd)

  cd.initialized = status
  if isready(cd.initUpChannel)
    tkst = take!(cd.initUpChannel)
    # @info "dumping stale cliq=$(cliq.index) status message $(tkst), replacing with $(status)"
  end
  put!(cd.initUpChannel, status)
  cond = getSolveCondition(cliq)
  notify(cond)
    # hack to avoid a race condition  -- remove with atomic lock logic upgrade
    sleep(0.1)
    notify(cond) # getSolveCondition(cliq)

  # TODO unlock
  unlockUpStatus!(cd)
  with_logger(logger) do
    tt = split(string(now()), 'T')[end]
    @info "$(tt) $(current_task()), cliq=$(cliq.index), notifyCliqUpInitStatus! -- unlocked, $(cd.initialized)"
  end

  nothing
end

function notifyCliqDownInitStatus!(cliq::TreeClique, status::Symbol; logger=ConsoleLogger())
  cdat = getCliqueData(cliq)
  with_logger(logger) do
    @info "$(now()) $(current_task()), cliq=$(cliq.index), notifyCliqDownInitStatus! -- pre-lock, new $(cdat.initialized)-->$(status)"
  end

  # take lock for atomic transaction
  lockDwnStatus!(cdat, cliq.index, logger=logger)

  cdat.initialized = status

  if isready(cdat.initDownChannel)
    content = take!(cdat.initDownChannel)
    with_logger(logger) do
      @info "dumping stale cliq=$(cliq.index) status message $(content), replacing with $(status)"
    end
  end
  put!(cdat.initDownChannel, status)
  notify(getSolveCondition(cliq))
    # hack to avoid a race condition
    sleep(0.1)
    notify(getSolveCondition(cliq))

  # unlock for others to proceed
  unlockDwnStatus!(cdat)
  with_logger(logger) do
    @info "$(now()), cliq=$(cliq.index), notifyCliqDownInitStatus! -- unlocked, $(getCliqStatus(cliq))"
  end

  # flush(logger.stream)

  nothing
end


"""
    $SIGNATURES

Return dictionary of all up belief messages currently in a Bayes `tree`.

Notes
- Returns `::Dict{Int,LikelihoodMessage}`

Related

getUpMsgs
"""
function getTreeCliqUpMsgsAll(tree::AbstractBayesTree)
  allUpMsgs = Dict{Int,LikelihoodMessage}()
  for (idx,cliq) in getCliques(tree)
    msgs = getUpMsgs(cliq)
    allUpMsgs[cliq.index] = LikelihoodMessage()
    for (lbl,msg) in msgs
      # TODO capture the inferred dimension as part of the upward propagation
      allUpMsgs[cliq.index].belief[lbl] = msg
    end
  end
  return allUpMsgs
end

"""
    $SIGNATURES

Convert tree up messages dictionary to a new dictionary relative to variables specific messages and their depth in the tree

Notes
- Return data in `TempUpMsgPlotting` format:
    Dict{Symbol,   -- is for variable label
     Vector{       -- multiple msgs for the same variable
      Symbol,      -- Clique index
      Int,         -- Depth in tree
      BTD          -- Belief estimate
      inferredDim  -- Information count
     }
"""
function stackCliqUpMsgsByVariable(tree::AbstractBayesTree,
                                   tmpmsgs::Dict{Int, LikelihoodMessage}  )
  #
  # start of the return data structure
  stack = TempUpMsgPlotting()

  # look at all the clique level data
  for (cidx,tmpmsg) in tmpmsgs
    # look at all variables up msg from each clique
    for (sym,msgdim) in tmpmsg.belief
      # create a new object for a particular variable if hasnt been seen before
      if !haskey(stack,sym)
        # FIXME this is an old message type
        stack[sym] = Vector{Tuple{Symbol, Int, BallTreeDensity, Float64}}()
      end
      # assemble metadata
      cliq = getCliques(tree,cidx)
      frt = getCliqFrontalVarIds(cliq)[1]
      # add this belief msg and meta data to vector of variable entry
      push!(stack[sym], (frt, getCliqDepth(tree, cliq),msgdim[1], msgdim[2]))
    end
  end

  return stack
end



## From other random places

"""
    $SIGNATURES

Return dictionary of down messages consisting of all frontal and separator beliefs of this clique.

Notes:
- Fetches numerical results from `subdfg` as dictated in `cliq`.
- return LikelihoodMessage
"""
function getCliqDownMsgsAfterDownSolve(subdfg::AbstractDFG, cliq::TreeClique)
  # Dict{Symbol, BallTreeDensity}
  # where the return msgs are contained
  container = LikelihoodMessage() # Dict{Symbol,BallTreeDensity}()

  # go through all msgs one by one
  for sym in getCliqAllVarIds(cliq)
    container.belief[sym] = TreeBelief( getVariable(subdfg, sym) )
  end

  # return the result
  return container
end



## =============================================================================
## DEPRECATED BELOW
## =============================================================================

#
