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
  # TODO remove once consolidation with upMsgs is done
  # putMsgUpInit!( cliq, LikelihoodMessage() )
  setCliqueMsgUp!(cda, LikelihoodMessage() )

  cda.dwnMsg = LikelihoodMessage()
  # cda.upInitMsgs = LikelihoodMessage()
  cda.downInitMsg = LikelihoodMessage()
  setCliqueStatus!(cliq, :null)
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
                        msgs::LikelihoodMessage;
                        tags::Vector{Symbol}=Symbol[])
  # add messages as priors to this sub factor graph
  msgfcts = DFGFactor[]
  svars = DFG.listVariables(subfg)
  tags__ = union(Symbol[:LIKELIHOODMESSAGE;], tags)
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
      fc = addFactor!(subfg, [msym], msgPrior, graphinit=false, tags=tags__)
      push!(msgfcts, fc)
    end
  end
  return msgfcts
end

function addMsgFactors!(subfg::AbstractDFG,
                        allmsgs::Dict{Int,LikelihoodMessage};
                        tags::Vector{Symbol}=Symbol[] )
  #
  allfcts = DFGFactor[]
  for (cliqid, msgs) in allmsgs
    # do each dict in array separately
    newfcts = addMsgFactors!(subfg, msgs, tags=tags)
    union!( allfcts, newfcts )
  end
  return allfcts
end


"""
    $SIGNATURES

Delete from the subgraph`::AbstractDFG` prior belief `msgs` that could/would be used
during clique inference.

DevNotes
- TODO make `::Vector{Symbol}` version.

Related

`addMsgFactors!`
"""
function deleteMsgFactors!(subfg::AbstractDFG,
                           fcts::Vector )
  #
  for fc in fcts
    deleteFactor!(subfg, fc.label)
  end
end
# deleteMsgFactors!(::LightDFG{SolverParams,DFGVariable,DFGFactor}, ::Array{DFGFactor{CommonConvWrapper{MsgPrior{BallTreeDensity}},1},1})



"""
    $SIGNATURES

Prepare the upward inference messages from clique to parent and return as `Dict{Symbol}`.

Notes
- Does not require tree message likelihood factors in subfg.
- Also see #579 regarding elimited likelihoods and priors.

DevNotes
- Consolidation in progress, part of #459
"""
function prepCliqInitMsgsUp(subfg::AbstractDFG,
                            cliq::TreeClique,
                            status::Symbol=getCliqueStatus(cliq);
                            logger=ConsoleLogger(),
                            duplicate::Bool=true )
  #
  # get the current clique status

  # construct init's up msg to place in parent from initialized separator variables
  msg = LikelihoodMessage(status)
  seps = getCliqSeparatorVarIds(cliq)
  with_logger(logger) do
    @info "prepCliqInitMsgsUp, seps=$seps"
  end
  for vid in seps
    var = DFG.getVariable(subfg, vid)
    var = duplicate ? deepcopy(var) : var
    if isInitialized(var)
      msg.belief[Symbol(var.label)] = TreeBelief(var)
    end
  end
  return msg
end




## =============================================================================
## Atomic messaging during init -- might be deprecated TODO
## =============================================================================




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

  setCliqueStatus!(cdat, status)

  putMsgDwnInitStatus!(cliq, status, logger)

  # unlock for others to proceed
  unlockDwnStatus!(cdat)
  with_logger(logger) do
    @info "$(now()), cliq=$(cliq.index), notifyCliqDownInitStatus! -- unlocked, $(getCliqueStatus(cliq))"
  end

  # flush(logger.stream)

  nothing
end




## =============================================================================
## Multimessage assemblies from multiple cliques
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
