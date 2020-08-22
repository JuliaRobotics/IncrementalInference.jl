

##==============================================================================
## LEGACY SUPPORT FOR ZMQ IN CAESAR
##==============================================================================

export listSolvekeys

export _evalType

# not sure if and where this is still being used
function _evalType(pt::String)::Type
    @error "_evalType has been deprecated, use DFG serialization methods instead."
    try
        getfield(Main, Symbol(pt))
    catch ex
        io = IOBuffer()
        showerror(io, ex, catch_backtrace())
        err = String(take!(io))
        error("_evalType: Unable to locate factor/distribution type '$pt' in main context (e.g. do a using on all your factor libraries). Please check that this factor type is loaded into main. Stack trace = $err")
    end
end


##==============================================================================
## Delete at end v0.15.x
##==============================================================================

# """
#     $SIGNATURES

# Get the latest down message from the parent node (without calculating anything).

# Notes
# - Different from down initialization messages that do calculate new values -- see `prepCliqInitMsgsDown!`.
# - Basically converts function `getDwnMsgs` from `Dict{Symbol,BallTreeDensity}` to `Dict{Symbol,Vector{BallTreeDensity}}`.

# Related

# getMsgInitDwnParent
# """
# function getMsgDwnParent(treel::AbstractBayesTree, cliq::TreeClique)
#   downmsgs = IntermediateMultiSiblingMessages()
#   prnts = getParent(treel, cliq)
#   if 0 < length(prnts)
#     prnt = prnts[1]
#     prntmsgs = getDwnMsgs(prnt)
#     for (key, bel) in prntmsgs
#       if !haskey(downmsgs, key)
#         downmsgs[key] = IntermediateSiblingMessages()
#       end
#       # TODO insert true inferred dim
#       push!(downmsgs[key], bel)
#     end
#   end

#   return downmsgs
# end


# function getMsgDwnThisInit(cdat::BayesTreeNodeData) 
#   iifdepwarn("#459 replace with getfetchCliqueMsgDown", :getMsgDwnThisInit)
#   return cdat.downInitMsg
# end

function addMsgFactors!(subfg::AbstractDFG,
                        msgs::Dict{Symbol, Vector{Tuple{BallTreeDensity, Float64}}} )
  # msgs::
  # add messages as priors to this sub factor graph
  # @warn "Tuple{KDE,Floa64} specific version of addMsgFactors! is deprecated, use LikelihoodMessage version instead."
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
