# init utils for tree based inference

export resetCliqSolve!


## =============================================================================
# short preamble funcions
## =============================================================================



convert(::Type{BallTreeDensity}, src::TreeBelief) = manikde!(src.val, src.bw[:,1], src.softtype)


## =============================================================================
# helper functions for tree message channels
## =============================================================================

"""
    $SIGNATURES

Reset the state of all variables in a clique to not initialized.

Notes
- resets numberical values to zeros.

Dev Notes
- TODO not all kde manifolds will initialize to zero.
- FIXME channels need to be consolidated
"""
function resetCliqSolve!(dfg::AbstractDFG,
                         treel::AbstractBayesTree,
                         cliq::TreeClique;
                         solveKey::Symbol=:default)
  #
  cda = getCliqueData(cliq)
  vars = getCliqVarIdsAll(cliq)
  for varis in vars
    resetVariable!(dfg, varis, solveKey=solveKey)
  end
  prnt = getParent(treel, cliq)
  if length(prnt) > 0
    putMsgUpInit!(prnt[1], cliq.index, LikelihoodMessage()) # TODO X putMsgUpInit!( cliq, cliq.index, LikelihoodMessage() )
  end
  cda.upMsg  = LikelihoodMessage()
  cda.dwnMsg = LikelihoodMessage()
  cda.upInitMsgs = Dict{Int, LikelihoodMessage}()
  cda.downInitMsg = LikelihoodMessage()
  setCliqStatus!(cliq, :null)
  setCliqDrawColor(cliq, "")
  return nothing
end

function resetCliqSolve!(dfg::AbstractDFG,
                         treel::AbstractBayesTree,
                         frt::Symbol;
                         solveKey::Symbol=:default  )
  #
  resetCliqSolve!(dfg, treel, getCliq(treel, frt), solveKey=solveKey)
end



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
- assume lower limit on number of particles is 5.
- messages from children stored in vector or dict.

Related

`deleteMsgFactors!`
"""
function addMsgFactors!(subfg::AbstractDFG,
                        msgs::LikelihoodMessage)
  # add messages as priors to this sub factor graph
  msgfcts = DFGFactor[]
  svars = DFG.listVariables(subfg)
  for (msym, belief_) in msgs.belief
    if msym in svars
      if 5 < size(belief_.val,2)
        # NOTE kde case assumes lower limit of 5 particles
        # manifold information is contained in the factor graph DFGVariable object
        kdePr = manikde!(belief_.val, belief_.bw[:,1], getManifolds(belief_.softtype))
        msgPrior = MsgPrior(kdePr, belief_.inferdim)
      elseif size(belief_.val, 2) == 1 && size(belief_.val, 1) == 1
        msgPrior =  MsgPrior(Normal(belief_.val[1], sqrt(belief_.bw[1])), belief_.inferdim)
      elseif size(belief_.val, 2) == 1 && 1 != size(belief_.val, 1)
        mvnorm = createMvNormal(belief_.val[:,1], belief_.bw)
        mvnorm != nothing ? nothing : (return DFGFactor[])
        msgPrior =  MsgPrior(mvnorm, belief_.inferdim)
      else
        error("Don't know what what to do with size(belief_.val)=$(size(belief_.val))")
      end
      fc = addFactor!(subfg, [msym], msgPrior, graphinit=false)
      push!(msgfcts, fc)
    end
  end
  return msgfcts
end

function addMsgFactors!(subfg::AbstractDFG,
                        allmsgs::Dict{Int,LikelihoodMessage} )
  #
  allfcts = DFGFactor[]
  for (cliqid, msgs) in allmsgs
    # do each dict in array separately
    newfcts = addMsgFactors!(subfg, msgs)
    union!( allfcts, newfcts )
  end
  return allfcts
end


"""
    $SIGNATURES

Delete from the subgraph`::AbstractDFG` prior belief `msgs` that could/would be used
during clique inference.

Related

`addMsgFactors!`
"""
function deleteMsgFactors!(subfg::AbstractDFG,
                           fcts::Vector{DFGFactor} )
  #
  for fc in fcts
    deleteFactor!(subfg, fc.label)
  end
end




## =============================================================================
## Atomic messaging during init -- might be deprecated TODO
## =============================================================================


"""
    $SIGNATURES

Update clique status and notify of the change

Notes
- Assumes users will lock the status state before getting status until after decision whether to update status.
- If so, only unlock after status and condition has been updated.

Dev Notes
- Should be made an atomic transaction
"""
function notifyCliqUpInitStatus!(cliq::TreeClique,
                                 status::Symbol;
                                 logger=ConsoleLogger() )
  #
  cd = getCliqueData(cliq)
  with_logger(logger) do
    tt = split(string(now()), 'T')[end]
    @info "$(tt) $(current_task()), cliq=$(cliq.index), notifyCliqUpInitStatus! -- pre-lock, $(cd.initialized)-->$(status)"
  end
  flush(logger.stream)

  # currently using a lock internally (hack message channels are consolidated)
  putMsgUpInitStatus!(cliq, status, logger)

  with_logger(logger) do
    tt = split(string(now()), 'T')[end]
    @info "$(tt) $(current_task()), cliq=$(cliq.index), notifyCliqUpInitStatus! -- unlocked, $(cd.initialized)"
  end

  nothing
end


function notifyCliqDownInitStatus!(cliq::TreeClique,
                                   status::Symbol;
                                   logger=ConsoleLogger() )
  #
  cdat = getCliqueData(cliq)
  with_logger(logger) do
    @info "$(now()) $(current_task()), cliq=$(cliq.index), notifyCliqDownInitStatus! -- pre-lock, new $(cdat.initialized)-->$(status)"
  end

  # take lock for atomic transaction
  lockDwnStatus!(cdat, cliq.index, logger=logger)

  cdat.initialized = status

  putMsgDwnInitStatus!(cliq, status, logger)

  # unlock for others to proceed
  unlockDwnStatus!(cdat)
  with_logger(logger) do
    @info "$(now()), cliq=$(cliq.index), notifyCliqDownInitStatus! -- unlocked, $(getCliqStatus(cliq))"
  end

  # flush(logger.stream)

  nothing
end




## =============================================================================
## Multimessage assemplies from multiple cliques
## =============================================================================


"""
    $SIGNATURES

Return dictionary of all up belief messages currently in a Bayes `tree`.

Notes
- Returns `::Dict{Int,LikelihoodMessage}`

Related

getUpMsgs
"""
function getTreeCliqUpMsgsAll(tree::AbstractBayesTree)
  allUpMsgs = Dict{Int,LikelihoodMessage}()
  for (idx,cliq) in getCliques(tree)
    msgs = getUpMsgs(cliq)
    allUpMsgs[cliq.index] = LikelihoodMessage()
    for (lbl,msg) in msgs
      # TODO capture the inferred dimension as part of the upward propagation
      allUpMsgs[cliq.index].belief[lbl] = msg
    end
  end
  return allUpMsgs
end

"""
    $SIGNATURES

Convert tree up messages dictionary to a new dictionary relative to variables specific messages and their depth in the tree

Notes
- Return data in `TempUpMsgPlotting` format:
    Dict{Symbol,   -- is for variable label
     Vector{       -- multiple msgs for the same variable
      Symbol,      -- Clique index
      Int,         -- Depth in tree
      BTD          -- Belief estimate
      inferredDim  -- Information count
     }
"""
function stackCliqUpMsgsByVariable(tree::AbstractBayesTree,
                                   tmpmsgs::Dict{Int, LikelihoodMessage}  )
  #
  # start of the return data structure
  stack = TempUpMsgPlotting()

  # look at all the clique level data
  for (cidx,tmpmsg) in tmpmsgs
    # look at all variables up msg from each clique
    for (sym,msgdim) in tmpmsg.belief
      # create a new object for a particular variable if hasnt been seen before
      if !haskey(stack,sym)
        # FIXME this is an old message type
        stack[sym] = Vector{Tuple{Symbol, Int, BallTreeDensity, Float64}}()
      end
      # assemble metadata
      cliq = getCliques(tree,cidx)
      frt = getCliqFrontalVarIds(cliq)[1]
      # add this belief msg and meta data to vector of variable entry
      push!(stack[sym], (frt, getCliqDepth(tree, cliq),msgdim[1], msgdim[2]))
    end
  end

  return stack
end



"""
    $SIGNATURES

Return dictionary of down messages consisting of all frontal and separator beliefs of this clique.

Notes:
- Fetches numerical results from `subdfg` as dictated in `cliq`.
- return LikelihoodMessage
"""
function getCliqDownMsgsAfterDownSolve(subdfg::AbstractDFG, cliq::TreeClique)
  # Dict{Symbol, BallTreeDensity}
  # where the return msgs are contained
  container = LikelihoodMessage() # Dict{Symbol,BallTreeDensity}()

  # go through all msgs one by one
  for sym in getCliqAllVarIds(cliq)
    container.belief[sym] = TreeBelief( getVariable(subdfg, sym) )
  end

  # return the result
  return container
end



#
