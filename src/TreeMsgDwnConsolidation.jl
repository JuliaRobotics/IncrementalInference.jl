
# these functions and their use in the code need to be consolidated to conclude #459

function putCliqueMsgDown!(cdata::BayesTreeNodeData, msg::LikelihoodMessage; from::Symbol=:nothing)
  @debug "putCliqueMsgDown! from=$(from)"
  cdata.dwnMsg = msg
end


function getfetchCliqueInitMsgDown(cdata::BayesTreeNodeData; from::Symbol=:nothing)
  @debug "getfetchCliqueInitMsgDown from=$(from)"
  return cdata.downInitMsg
end
@deprecate getfetchCliqueMsgDown(cdata::BayesTreeNodeData; from::Symbol=:nothing) getfetchCliqueInitMsgDown(cdata, from=from)

function putCliqueInitMsgDown!(cdata::BayesTreeNodeData, initmsg::LikelihoodMessage)
  cdata.downInitMsg = initmsg
  nothing
end


"""
    $(SIGNATURES)

Set the downward passing message for Bayes (Junction) tree clique `cliql`.
"""
function putMsgDwnThis!(cliql::TreeClique, msgs::LikelihoodMessage)
  iifdepwarn("#459 replace with putCliqueMsgDown!", :putMsgDwnThis!)
  getCliqueData(cliql).dwnMsg = msgs
end
putMsgDwnThis!(csmc::CliqStateMachineContainer, msgs::LikelihoodMessage) = putMsgDwnThis!(csmc.cliq, msgs)  # NOTE, old, csmc.msgsDown = msgs


"""
    $(SIGNATURES)

Return the last down message stored in `cliq` of Bayes (Junction) tree.
"""
fetchMsgDwnThis(cliql::TreeClique) = getCliqueData(cliql).dwnMsg
fetchMsgDwnThis(csmc::CliqStateMachineContainer) = fetchMsgDwnThis(csmc.cliq)
fetchMsgDwnThis(btl::AbstractBayesTree, sym::Symbol) = fetchMsgDwnThis(getClique(btl, sym))

getMsgDwnInitChannel_(cdat::BayesTreeNodeData) = cdat.initDownChannel

getMsgDwnThisInit(cliq::TreeClique) = getMsgDwnThisInit(getCliqueData(cliq))

getMsgDwnInitChannel_(cliq::TreeClique) = getMsgDwnInitChannel_(getCliqueData(cliq))
fetchMsgDwnInit(cliq::TreeClique) = fetch(getMsgDwnInitChannel_(cliq))


function blockMsgDwnUntilStatus(cliq::TreeClique, status::CliqStatus)
  while fetchMsgDwnInit(cliq).status != status
    wait(getSolveCondition(cliq))
  end
  nothing
end


# FIXME will be consolidated as part of 459
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

Calculate a new down message from the parent.
"""
function getMsgInitDwnParent(treel::AbstractBayesTree, cliq::TreeClique; logger=SimpleLogger(stdout))
  # FIXME drop IntermediateMultiSiblingMessages and use only LikelihoodMessage
  # check if any msgs should be multiplied together for the same variable
  # msgspervar = LikelihoodMessage() # or maybe Dict{Int, LikelihoodMessage}()

  # get the parent
  prnt = getParent(treel, cliq)[1]
  msgspervar = IntermediateMultiSiblingMessages()

  # get the current messages ~~stored in~~ [going to] the parent (pull model #674)
  prntmsgs = getMsgsUpInitChildren(treel, prnt, TreeBelief, skip=[cliq.index;])     # FIXME, post #459 calls?

  with_logger(logger) do
    @info "prnt $(prnt.index), getMsgInitDwnParent -- msg ids::Int=$(collect(keys(prntmsgs)))"
  end

  for (msgcliqid, msgs) in prntmsgs
    with_logger(logger) do
      @info "getMsgInitDwnParent -- msgcliqid=$msgcliqid, msgs.belief=$(collect(keys(msgs.belief)))"
    end
    for (msgsym, msg) in msgs.belief
      if !haskey(msgspervar, msgsym)
        # there will be an entire list...
        msgspervar[msgsym] = IntermediateSiblingMessages()
      end
      with_logger(logger) do
        @info "getMsgInitDwnParent -- msgcliqid=$(msgcliqid), msgsym $(msgsym), inferdim=$(msg.inferdim)"
      end
      push!(msgspervar[msgsym], msg)
    end
  end

  return msgspervar
end



## =============================================================
