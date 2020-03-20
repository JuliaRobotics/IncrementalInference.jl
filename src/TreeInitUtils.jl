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
                        msgs::Dict{Symbol, Vector{Tuple{BallTreeDensity, Float64}}} )::Vector{DFGFactor}
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
                        allmsgs::Dict{Int,LikelihoodMessage} )::Vector{DFGFactor}
  #
  allfcts = DFGFactor[]
  for (cliqid, msgs) in allmsgs
    # do each dict in array separately
    newfcts = addMsgFactors!(subfg, msgs)
    union!( allfcts, newfcts )
  end
  return allfcts
end


# Consolidate with nonparametric addMsgFactors! ?
function addMsgFactors_Parametric!(subfg::AbstractDFG,
                                   msgs::BeliefMessage)::Vector{DFGFactor}
  # add messages as priors to this sub factor graph
  msgfcts = DFGFactor[]
  svars = DFG.listVariables(subfg)
  for (msym, belief) = (msgs.belief)
    if msym in svars
      #TODO covaraince
      #TODO Maybe always use MvNormal
      if size(belief.val)[1] == 1
        msgPrior =  MsgPrior(Normal(belief.val[1], sqrt(belief.bw[1])), belief.inferdim)
      else
        #FIXME a hack to make matrix Hermitian
        covar = Symmetric(belief.bw + 1e-5I)
        try
          MvNormal(belief.val[:,1], covar)
        catch er
          @error er "MvNormal Failed with:" covar
          return DFGFactor[]
        end
        msgPrior =  MsgPrior(MvNormal(belief.val[:,1], covar), belief.inferdim)
      end
      fc = addFactor!(subfg, [msym], msgPrior, graphinit=false)
      push!(msgfcts, fc)
    end
  end
  return msgfcts
end



"""
    $SIGNATURES

Delete from the subgraph`::AbstractDFG` prior belief `msgs` that could/would be used
during clique inference.

Related

`addMsgFactors!`
"""
function deleteMsgFactors!(subfg::AbstractDFG,
                           fcts::Vector{DFGFactor})
  #
  for fc in fcts
    deleteFactor!(subfg, fc.label)
  end
end




## =============================================================================
# tree belief message accessors for nonparameteric in this section
## =============================================================================


"""
    $SIGNATURES

Based on a push model from child cliques that should have already completed their computation.
"""
getCliqInitUpMsgs(cliq::TreeClique)::Dict{Int, LikelihoodMessage} = getCliqueData(cliq).upInitMsgs


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
  getCliqueData(cliql).upMsg = msgs
  nothing
end

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
function setDwnMsg!(cliql::TreeClique, msgs::LikelihoodMessage) #Dict{Symbol, BallTreeDensity}
  getCliqueData(cliql).dwnMsg = msgs
end

"""
    $(SIGNATURES)

Return the last down message stored in `cliq` of Bayes (Junction) tree.
"""
getDwnMsgs(cliql::TreeClique) = getCliqueData(cliql).dwnMsg
getDwnMsgs(btl::AbstractBayesTree, sym::Symbol) = getDwnMsgs(getCliq(btl, sym))






## =============================================================================
## DEPRECATED BELOW
## =============================================================================


# """
#     $(SIGNATURES)
#
# Return the last down message stored in `cliq` of Bayes (Junction) tree.
# """
# getCliqMsgsDown(cliql::TreeClique) = getDwnMsgs(cliql)

# """
#     $(SIGNATURES)
#
# Return the last up message stored in `cliq` of Bayes (Junction) tree.
# """
# getCliqMsgsUp(cliql::TreeClique) = upMsg(cliql)
# getCliqMsgsUp(treel::AbstractBayesTree, frt::Symbol) = getCliqMsgsUp(getCliq(treel, frt))


# function setUpMsg!(cliql::TreeClique, msgs::Dict{Symbol, BallTreeDensity})
#   @error "setUpMsg!, use inferred dimension version instead"
#   getCliqueData(cliql).upMsg = msgs
# end
#
