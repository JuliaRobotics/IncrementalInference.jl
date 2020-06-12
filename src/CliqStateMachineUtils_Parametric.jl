#TODO move to TreeMessageUtils.jl after development and merge
# retun ::Vector{DFGFactor}
function addMsgFactors_Parametric!(subfg::AbstractDFG,
                                   msgs::LikelihoodMessage)
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
