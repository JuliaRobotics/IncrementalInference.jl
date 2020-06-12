#TODO move to TreeMessageUtils.jl after development and merge
function addMsgFactors_Parametric!(subfg::AbstractDFG,
                                   msgs::LikelihoodMessage)::Vector{DFGFactor}
  # add messages as priors to this sub factor graph
  msgfcts = DFGFactor[]
  svars = DFG.listVariables(subfg)
# <<<<<<< Updated upstream
  for (msym, belief) = (msgs.belief)
    if msym in svars
      #TODO covaraince
      #TODO Maybe always use MvNormal
      if size(belief.val)[1] == 1
        msgPrior =  MsgPrior(Normal(belief.val[1], sqrt(belief.bw[1])), belief.inferdim)
      else
        mvnorm = createMvNormal(belief.val[:,1], belief.bw)
        mvnorm == nothing &&
          (return DFGFactor[])
        msgPrior =  MsgPrior(mvnorm, belief.inferdim)
      end
      fc = addFactor!(subfg, [msym], msgPrior, graphinit=false)
      push!(msgfcts, fc)
    end
# =======
#   #convert to Normal or MvNormal
#   varOrder = msgs.cobelief.varlbl
#   Σ = msgs.cobelief.Σ
#   μ = msgs.cobelief.μ
#   if length(μ) == 1
#     dist = Normal(μ[1], sqrt(Σ[1]))
#   else
#     dist = MvNormal(μ, Σ)
#   end
#
#   #TODO add type of factor to message, maybe even send constucted function
#   #TODO inferredDim
#   if length(varOrder) == 1
#     @info "Adding belief msg prior with $dist on $varOrder"
#     fc = addFactor!(subfg, varOrder, MsgPrior(dist, 0.0), graphinit=false)
#   elseif length(varOrder) == 2
#     @error "this only works for linear conditional"
#     fc = addFactor!(subfg, varOrder, LinearConditional(dist), graphinit=false)
#   else
#     error("Oops, not supported")
#   end
#   push!(msgfcts, fc)
#
#   # for (msym, belief) = (msgs.belief)
#   #   if msym in svars
#   #     #TODO covaraince
#   #     #TODO Maybe always use MvNormal
#   #     if size(belief.val)[1] == 1
#   #       msgPrior =  MsgPrior(Normal(belief.val[1], belief.bw[1]), belief.inferdim)
#   #     else
#   #       msgPrior =  MsgPrior(MvNormal(belief.val[:,1], belief.bw), belief.inferdim)
#   #     end
#   #     fc = addFactor!(subfg, [msym], msgPrior, graphinit=false)
#   #     push!(msgfcts, fc)
#   #   end
# >>>>>>> Stashed changes
  end
  return msgfcts
end

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
