
export
  getCliqueStatus,
  setCliqueStatus!,
  getSolveCondition

# Reguler accessors
export
  getMsgUpThis,
  getMsgsUpChildren
  # putMsgUpThis!

export
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
## Regular up and down Message Registers/Channels, getters and setters
## =============================================================================

## =============================================================================
## Consolidating INIT up and down Message Registers/Channels, getters and setters
## =============================================================================

# FIXME DEPRECATE to use only MsgUpThis register
getMsgUpChannel(tree::BayesTree, edge) = tree.messages[edge.index].upMsg
getMsgUpChannel(tree::MetaBayesTree, edge) = MetaGraphs.get_prop(tree.bt, edge, :upMsg)

# Consolidate to only use this
getMsgUpChannel(cdat::BayesTreeNodeData) = cdat.upMsgChannel
getMsgUpChannel(cliq::TreeClique) = getMsgUpChannel(getCliqueData(cliq))




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

Put a belief message on the up tree message channel `edge`. Blocks until a take! is performed by a different task.
"""
function putBeliefMessageUp!(tree::AbstractBayesTree, edge, beliefMsg::LikelihoodMessage)
  # Blocks until data is available.
  put!(getMsgUpChannel(tree, edge), beliefMsg)
  return beliefMsg
end


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
function prepPutCliqueStatusMsgUp!( csmc::CliqStateMachineContainer,
                                    status::Symbol=getCliqueStatus(csmc.cliq);
                                    dfg::AbstractDFG=csmc.cliqSubFg,
                                    upmsg=prepCliqueMsgUpConsolidated(dfg, csmc.cliq, status)  )
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
function getMsgsUpChildren( treel::AbstractBayesTree,
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
