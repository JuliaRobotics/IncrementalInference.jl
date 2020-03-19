# init utils for tree based inference




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


## MUST CONSOLIDATE
function addMsgFactors_Parametric!(subfg::G,
                                   msgs::BeliefMessage)::Vector{DFGFactor} where G <: AbstractDFG
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




### DEPRECATED BELOW




# function addMsgFactors!(subfg::G, msgs::LikelihoodMessage) where G <: AbstractDFG
#   @warn "addMsgFactors! use LikelihoodMessage format instead"
#
#   # assemble new temporary list
#   tmpmsgs = LikelihoodMessage()
#   for (id, val) in msgs
#     tmpmsgs.belief[id] = (val, 0.0)
#   end
#
#   # call intended function with temporary message work around
#   return addMsgFactors!(subfg, tmpmsgs)
# end
