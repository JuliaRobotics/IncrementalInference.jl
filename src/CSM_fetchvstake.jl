# Must consolidate as part of #855 , fetch vs take!

# Regular accessors
export
  fetchMsgUpThis,
  fetchMsgsUpChildren,
  fetchChildrenStatusUp,
  getMsgsUpInitChildren




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
    $(SIGNATURES)

Return the last up message stored in This `cliq` of the Bayes (Junction) tree.
"""
fetchMsgUpThis(cdat::BayesTreeNodeData) = fetch(getMsgUpChannel(cdat))  # cdat.upMsg    # TODO rename to fetchMsgUp
fetchMsgUpThis(cliql::TreeClique) = fetchMsgUpThis(getCliqueData(cliql))
fetchMsgUpThis(btl::AbstractBayesTree, frontal::Symbol) = fetchMsgUpThis(getClique(btl, frontal))



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
- FIXME update refactor to fetch or take, #855
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


function fetchMsgsUpChildren(csmc::CliqStateMachineContainer,
                            ::Type{TreeBelief}=TreeBelief )
  #
  # TODO, replace with single channel stored in csmcs or cliques
  fetchMsgsUpChildren(csmc.tree, csmc.cliq, TreeBelief)
end

## ====================================================================================
## TODO Deprecate below
## ====================================================================================


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


