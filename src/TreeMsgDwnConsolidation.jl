
# these functions and their use in the code need to be consolidated to conclude #459

# TODO consolidate
export
  fetchDwnMsgsThis,
  getMsgDwnThisInit,
  getMsgDownParent,
  getMsgDwnInitChannel_


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



## ============================================================================
## .initDownChannel, MUST BE REMOVED
## ============================================================================


@deprecate putMsgDwnInitChannel!(btnd::BayesTreeNodeData, msg::LikelihoodMessage) putDwnMsgConsolidated!(btnd, msg)
@deprecate getMsgDwnInitChannel_(btnd::BayesTreeNodeData) getDwnMsgConsolidated(btnd)
# getMsgDwnInitChannel_(btnd::BayesTreeNodeData) = btnd.initDownChannel

getMsgDwnInitChannel_(cliq::TreeClique) = getMsgDwnInitChannel_(getCliqueData(cliq))
fetchMsgDwnInit(cliq::TreeClique) = fetch(getMsgDwnInitChannel_(cliq))



# FIXME OLD must be consolidated as part of 459
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

## ============================================================================
## .downInitMsg, MUST BE REMOVED
## ============================================================================


function getfetchCliqueInitMsgDown(cdata::BayesTreeNodeData; from::Symbol=:nothing)
  @debug "getfetchCliqueInitMsgDown from=$(from)"
  return cdata.downInitMsg
end
# getMsgDwnThisInit(cliq::TreeClique) = getMsgDwnThisInit(getCliqueData(cliq)) # WHAT ???

function putCliqueInitMsgDown!(cdata::BayesTreeNodeData, initmsg::LikelihoodMessage)
  cdata.downInitMsg = initmsg
  nothing
end



## ============================================================================
## .dwnMsg, MUST BE REMOVED
## ============================================================================

"""
    $(SIGNATURES)

Return the last down message stored in `cliq` of Bayes (Junction) tree.
"""
fetchMsgDwnThis(cliql::TreeClique) = getCliqueData(cliql).dwnMsg
fetchMsgDwnThis(csmc::CliqStateMachineContainer) = fetchMsgDwnThis(csmc.cliq)
fetchMsgDwnThis(btl::AbstractBayesTree, sym::Symbol) = fetchMsgDwnThis(getClique(btl, sym))



"""
$(SIGNATURES)

Set the downward passing message for Bayes (Junction) tree clique `cliql`.
"""  
function putMsgDwnThis!(cdata::BayesTreeNodeData, msg::LikelihoodMessage; from::Symbol=:nothing)
  @debug "putMsgDwnThis! from=$(from)"
  cdata.dwnMsg = msg
end  
function putMsgDwnThis!(cliql::TreeClique, msgs::LikelihoodMessage)
  getCliqueData(cliql).dwnMsg = msgs
end  
putMsgDwnThis!(csmc::CliqStateMachineContainer, msgs::LikelihoodMessage) = putMsgDwnThis!(csmc.cliq, msgs)  # NOTE, old, csmc.msgsDown = msgs
  


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



function iifdepwarn(msg, funcsym; maxlog=nothing)
  @logmsg(
      Base.CoreLogging.Warn,
      msg,
      _module=begin
          bt = backtrace()
          frame, caller = Base.firstcaller(bt, funcsym)
          # TODO: Is it reasonable to attribute callers without linfo to Core?
          caller.linfo isa Core.MethodInstance ? caller.linfo.def.module : Core
      end,
      _file=String(caller.file),
      _line=caller.line,
      _id=(frame,funcsym),
      _group=:iifdepwarn,
      caller=caller,
      short_stacktrace=stacktrace(bt)[7:9],
      maxlog=maxlog
  )
  nothing
end

function Base.getproperty(obj::BayesTreeNodeData, sym::Symbol)
  if sym == :dwnMsg
    # iifdepwarn("#459 get dwnMsg", :getproperty)
  elseif sym == :downInitMsg
    # iifdepwarn("#459 get downInitMsg", :getproperty)
  elseif sym == :initDownChannel
    # iifdepwarn("#459 get initDownChannel", :getproperty)
  end
  return getfield(obj, sym)
end

function Base.setproperty!(obj::BayesTreeNodeData, sym::Symbol, val)
  if sym == :dwnMsg
    # iifdepwarn("#459 set dwnMsg", :setproperty!)
  elseif sym == :downInitMsg
    # iifdepwarn("#459 set downInitMsg", :setproperty!)
  elseif sym == :initDownChannel
    # iifdepwarn("#459 set initDownChannel", :setproperty!)
  end
  return setfield!(obj, sym, convert(fieldtype(typeof(obj), sym), val))
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
  lockDwnStatus!(cdat, cliq.index, logger=logger)

  setCliqueStatus!(cdat, status)

  putMsgDwnInitStatus!(cliq, status, logger)

  # unlock for others to proceed
  unlockDwnStatus!(cdat)
    with_logger(logger) do
    @info "$(now()), cliq=$(cliq.index), notifyCliqDownInitStatus! -- unlocked, $(getCliqueStatus(cliq))"
  end

  # flush(logger.stream)

  nothing
end





## =============================================================
