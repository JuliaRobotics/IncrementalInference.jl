

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
## Delete at end v0.16.x
##==============================================================================


"""
    $SIGNATURES

Set the marginalized status of a clique.
"""
function setCliqAsMarginalized!(cliq::TreeClique, status::Bool)
  @warn "Busy deprecating setCliqAsMarginalized!, use CSM features directly"
  if status
    setCliqueStatus!(cliq, :marginalized)
  else
    if getCliqueData(cliq).initialized == :marginalized
      @info "Reverting clique $(cliq.index) to assumed :downsolved status"
      setCliqueStatus!(cliq, :downsolved)
    else
      error("Unknown clique de-marginalization requist for clique $(cliq.index), current status: $(cliq.initialized)")
    end
  end
end

"""
    $SIGNATURES

Run through entire tree and set cliques as marginalized if all clique variables are marginalized.

Notes:
- TODO can be made fully parallel, consider converting for use with `@threads` `for`.
"""
function updateTreeCliquesAsMarginalizedFromVars!(fgl::AbstractDFG, tree::AbstractBayesTree)::Nothing
  @warn "Busy deprecating updateTreeCliquesAsMarginalizedFromVars!, use CSM features directly"
  for (clid, cliq) in getCliques(tree)
    if isCliqMarginalizedFromVars(fgl, cliq)
      setCliqAsMarginalized!(cliq, true)
    end
  end
  nothing
end


##==============================================================================
## Delete at end v0.15.x
##==============================================================================


@deprecate putCliqueMsgDown!(cdata::BayesTreeNodeData, msg::LikelihoodMessage; from::Symbol=:nothing) putMsgDwnThis!(cdata, msg, from=from)

# """
#     $SIGNATURES

# Initialization downward message passing is different from regular inference since
# it is possible that none of the child cliq variables have been initialized.

# Notes
# - init upward msgs are individually stored in child cliques (pull model good).
# - fresh product of overlapping beliefs are calculated on each function call.
# - Assumed that `prnt` of siblings

# Dev Notes
# - This should be the initialization cycle of parent, build up bit by bit...
# """
# function prepCliqInitMsgsDown!(fgl::AbstractDFG,
#                                tree::AbstractBayesTree,
#                                prnt::TreeClique,
#                                cliq::TreeClique;
#                                logger=ConsoleLogger() )
#   #
#   # tt = split(string(now()), 'T')[end]
#   # with_logger(logger) do
#   #   @info "$(tt) prnt $(prnt.index), prepCliqInitMsgsDown! -- with cliq $(cliq.index)"
#   # end
  
#   # FIXME use only LikelihoodMessage
#   # check if any msgs should be multiplied together for the same variable
#   # msgspervar = LikelihoodMessage() # or maybe Dict{Int, LikelihoodMessage}()
#   msgspervar = getMsgInitDwnParent(tree, cliq, logger=logger)
#   # reference to default dict location
#   #JT 459 products = getMsgDwnThisInit(prnt)
#   products = getfetchCliqueInitMsgDown(prnt.data, from=:getMsgDwnThisInit) |> deepcopy
  
#   ## TODO use parent factors too
#   # intersect with the asking clique's separator variables
#   condenseDownMsgsProductPrntFactors!(fgl, products, msgspervar, prnt, cliq, logger)
  
#   # with_logger(logger) do
#   #   @info "cliq $(prnt.index), prepCliqInitMsgsDown! -- vars fw/ down msgs=$(collect(keys(msgspervar)))"
#   # end
#   # flush(logger.stream)

#   # remove msgs that have no data
#   rmlist = Symbol[]
#   for (prsym,belmsg) in products.belief
#     if belmsg.inferdim < 1e-10
#       # no information so remove
#       push!(rmlist, prsym)
#     end
#   end
#   with_logger(logger) do
#     @info "cliq $(prnt.index), prepCliqInitMsgsDown! -- rmlist, no inferdim, keys=$(rmlist)"
#   end
#   for pr in rmlist
#     delete!(products.belief, pr)
#   end

#   with_logger(logger) do
#     @info "cliq $(prnt.index), prepCliqInitMsgsDown! -- product keys=$(collect(keys(products.belief)))"
#   end

#   # now put the newly computed message in the appropriate container
#   # FIXME THIS IS A PUSH MODEL, see #674 -- make pull model
#   putCliqueInitMsgDown!(getCliqueData(prnt), products)

#   return products
# end


# used during nonparametric CK preparation, when information from multiple siblings must be shared together
# const IntermediateSiblingMessagesTB{T} = Vector{TreeBelief{T}}
# const IntermediateMultiSiblingMessagesTB{T} = Dict{Symbol, Vector{TreeBelief{T}}}


# FIXME, better standardize intermediate types
# can be replaced by Vector{TreeBelief}
# const IntermediateSiblingMessages = Vector{Tuple{BallTreeDensity,Float64}}
# const IntermediateMultiSiblingMessages = Dict{Symbol, IntermediateSiblingMessages}


# Helper function for prepCliqInitMsgsDown!
# populate products with products of upward messages
function condenseDownMsgsProductOnly!(fgl::AbstractDFG,
                                      products::LikelihoodMessage,
                                      msgspervar::Dict{Symbol, <:AbstractVector}  )
  #
  error("condenseDownMsgsProductOnly!(::AbstractDFG,::LikelihoodMessage, ::IntermediateMultiSiblingMessages) is obsolete")
  # multiply multiple messages together
  for (msgsym, msgsBo) in msgspervar
    # check if this particular down message requires msgsym
    if exists(fgl, msgsym) # DFG.hasVariable(fgl, msgsym)
      if length(msgspervar[msgsym]) > 1
        msgs = getindex.(msgsBo, 1)
        haspars = 0.0
        for mb in msgsBo, val in mb[2]
          haspars += val
        end
        products[msgsym] = (manifoldProduct(msgs, getManifolds(fgl, msgsym)), haspars)
      else
        # transfer if only have a single belief
        products[msgsym] = (msgsBo[1][1], msgsBo[1][2])
      end
    else
      # not required, therefore remove from message to avoid confusion
      if haskey(products, msgsym)
        delete!(products, msgsym)
      end
    end
  end
  nothing
end

@deprecate getfetchCliqueMsgDown(cdata::BayesTreeNodeData; from::Symbol=:nothing) getfetchCliqueInitMsgDown(cdata, from=from)

export getIdx

"""
    $SIGNATURES

Return interger index of desired variable element.

Example
-------
```julia
pp = RoME.Point2()
getIdx(pp, :posY) # = 2
```

Internal Notes
--------------
- uses number i < 100 for index number, and
- uses +100 offsets to track the minibatch number of the requested dimension
"""
function getIdx(pp::Tuple,
                sym::Symbol,
                i::Int=0)
  #
  error("getIdx is obsolete, use DistributedFactorGraphs objects/methods instead.")
  i-=100
  for p in pp
    i,j = getIdx(p, sym, i)
    if i > 0
      return i, j
    end
  end
  return i,-1
end
getIdx(pp::Symbol, sym::Symbol, i::Int=0) = pp==sym ? (abs(i)%100+1, div(abs(i)-100,100)) : (i-1, div(abs(i)-100,100))
function getIdx(pp::InferenceVariable, sym::Symbol, i::Int=0)
  return getIdx(pp.dimtype, sym)
end


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

@deprecate fetchDataElement(dfg::AbstractDFG, varsym::Symbol, lbl::Symbol) fetchDataJSON(dfg, varsym, lbl)


function addMsgFactors!(subfg::AbstractDFG,
                        msgs::Dict{Symbol, <:AbstractVector} )  #Dict{Symbol, Vector{Tuple{BallTreeDensity, Float64}}} )
  # msgs::
  # add messages as priors to this sub factor graph
  @warn "Tuple{KDE,Floa64} specific version of addMsgFactors! is deprecated, use LikelihoodMessage version instead."
  msgfcts = DFGFactor[]
  svars = DFG.listVariables(subfg)
  for (msym, dms) in msgs
    for treebelief in dms
      if msym in svars
        # TODO should be on manifold prior, not just generic euclidean prior -- okay since variable on manifold, but not for long term
        fc = addFactor!(subfg, [msym], MsgPrior(manikde!(treebelief), treebelief.inferdim), graphinit=false)
        push!(msgfcts, fc)
      end
    end
  end
  return msgfcts
end
