"""
    $SIGNATURES

EXPERIMENTAL: Init and start state machine for parametric solve.
"""
function initStartCliqStateMachineParametric!(dfg::G,
                                              tree::AbstractBayesTree,
                                              cliq::TreeClique,
                                              cliqKey::Int; # moved to end during #459
                                              oldcliqdata::BayesTreeNodeData=BayesTreeNodeData(),
                                              verbose::Bool=false,
                                              drawtree::Bool=false,
                                              show::Bool=false,
                                              incremental::Bool=true,
                                              limititers::Int=-1,
                                              upsolve::Bool=true,
                                              downsolve::Bool=true,
                                              recordhistory::Bool=false,
                                              delay::Bool=false,
                                              logger::SimpleLogger=SimpleLogger(Base.stdout)) where {G <: AbstractDFG, AL <: AbstractLogger}

  # TODO Taking out, use tree and messages for operations involving children and parents
  children = TreeClique[]#getChildren(tree, cliq)
  prnt = TreeClique[]#getParent(tree, cliq)

  destType = (G <: InMemoryDFGTypes) ? G : InMemDFGType

  csmc = CliqStateMachineContainer(dfg, initfg(destType, solverParams=getSolverParams(dfg)),
                                    tree, cliq, prnt, children,
                                    incremental, drawtree, downsolve, delay,
                                    getSolverParams(dfg), Dict{Symbol,String}(), oldcliqdata, logger,
                                    cliqKey )
  #
  # nxt = upsolve ? testCliqCanRecycled_ParametricStateMachine : (downsolve ? testCliqCanRecycled_ParametricStateMachine : error("must attempt either up or down solve"))
  nxt = buildCliqSubgraph_ParametricStateMachine

  statemachine = StateMachine{CliqStateMachineContainer}(next=nxt)


  # store statemachine and csmc in task
  if dfg.solverParams.dbg || recordhistory
    task_local_storage(:statemachine, statemachine)
    task_local_storage(:csmc, csmc)
  end

  while statemachine(csmc, verbose=verbose, iterlimit=limititers, recordhistory=recordhistory); end
  statemachine.history
end

"""
    $SIGNATURES

Build a sub factor graph for clique variables from the larger factor graph.

Notes
- Parametric State machine function nr.1
"""
function buildCliqSubgraph_ParametricStateMachine(csmc::CliqStateMachineContainer)
  # build a local subgraph for inference operations
  syms = getCliqAllVarIds(csmc.cliq)

  infocsm(csmc, "Par-1, build subgraph syms=$(syms)")

  frontsyms = getCliqFrontalVarIds(csmc.cliq)
  sepsyms = getCliqSeparatorVarIds(csmc.cliq)
  buildCliqSubgraph!(csmc.cliqSubFg, csmc.dfg, frontsyms, sepsyms)

  #TODO remove, just as a sanity check if any priors remains on seperators
  removedIds = removeSeparatorPriorsFromSubgraph!(csmc.cliqSubFg, csmc.cliq)
  length(removedIds) > 0 && @error "removeSeparatorPriorsFromSubgraph, removed priors should not happen"


  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_build")

  # go to 2 wait for up
  return waitForUp_ParametricStateMachine
end

"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 2
- `take!` model #855
"""
function waitForUp_ParametricStateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "Par-2, wait for up messages of needed")

  setCliqDrawColor(csmc.cliq, "purple") #TODO don't know if this is correct color

  # JT empty upRx buffer to save messages, TODO It may be ok not to empty 
  beliefMessages = empty!(messages(csmc.cliq).upRx)

  # fetch messages from edges
  @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
    @async begin
      thisEdge = isa(e,Graphs.Edge) ? e.index : e
      @info "$(csmc.cliq.index): take! on edge $thisEdge"
      # Blocks until data is available. -- take! model
      beliefMsg = takeBeliefMessageUp!(csmc.tree, e)
      beliefMessages[thisEdge] = beliefMsg
      @info "$(csmc.cliq.index): Belief message received with status $(beliefMsg.status)"
    end
  end

  # get all statuses from messages
  all_child_status = map(msg -> msg.status, values(beliefMessages))

  # Main Branching happens here - all up messages received

  # If one up error is received propagate ERROR_STATUS 
  if :ERROR_STATUS in all_child_status

    setCliqDrawColor(csmc.cliq, "red")

    for e in getEdgesParent(csmc.tree, csmc.cliq)
      @info "Par-2, $(csmc.cliq.index): propagate up error on edge $(isa(e,Graphs.Edge) ? e.index : e)"
      putBeliefMessageUp!(csmc.tree, e, LikelihoodMessage(status=:ERROR_STATUS))
    end
    #if its a root, propagate error down
    #FIXME rather check if no parents with function (hasParents or isRoot)
    if length(getParent(csmc.tree, csmc.cliq)) == 0
      @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
        @info "Par-2 Root $(csmc.cliq.index): propagate down error on edge $(isa(e,Graphs.Edge) ? e.index : e)"
        @async putBeliefMessageDown!(csmc.tree, e, LikelihoodMessage(status=:ERROR_STATUS))
      end
      @error "Par-2 $(csmc.cliq.index): Exit with error state"
      return IncrementalInference.exitStateMachine
    end
    return waitForDown_X_StateMachine

  # just guessing some other states here as an example
  # elseif :UP_INITIALIZED in all_child_status
  #   return ...
  
  elseif all(all_child_status .== :UPSOLVED)
    return solveUp_ParametricStateMachine

  else
    error("waitForUp State Error: Unknown transision.")
  end
  
end

# Graph.jl does not have an in_edges function for a GenericIncidenceList, so extending here.
function Graphs.in_edges(vert::V, gr::GenericIncidenceList{V, Edge{V}, Vector{V}}) where {V}
  inclist = gr.inclist
  targid = vert.index
  inlist = Edge{V}[]
  for edgelist in inclist
    for ed in edgelist
      if ed.target.index == targid
        push!(inlist, ed)
      end
    end
  end
  return inlist
end

"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 3
"""
function solveUp_ParametricStateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "Par-3, Solving Up")

  setCliqDrawColor(csmc.cliq, "red")
  # csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  #TODO maybe change to symbols
  msgfcts = DFGFactor[]
  # LITTLE WEIRD get previously set up msgs (stored in this clique)
    # FIXME, fetch message buffered in channels
  # see #855
  for (idx,upmsg) in messages(csmc.cliq).upRx #get cached messages taken from children saved in this clique
    #TODO remove temp msgfcts container
    append!(msgfcts, addMsgFactors!(csmc.cliqSubFg, upmsg, UpwardPass) ) # addMsgFactors_Parametric!
  end
  @info "length mgsfcts=$(length(msgfcts))"
  infocsm(csmc, "length mgsfcts=$(length(msgfcts))")

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_beforeupsolve")

  vardict, result, varIds, Σ = solveFactorGraphParametric(csmc.cliqSubFg)

  @info "$(csmc.cliq.index) vars $(keys(varIds))"
  # @info "$(csmc.cliq.index) Σ $(Σ)"
  # Pack all results in variables
  # FIXME test f_converged, ls_success, confirm convergence check
  if result.f_converged || result.g_converged
    @info "$(csmc.cliq.index): subfg optim converged updating variables"
    for (v,val) in vardict
      vnd = getSolverData(getVariable(csmc.cliqSubFg, v), :parametric)
      # fill in the variable node data value
      @info "$(csmc.cliq.index) up: updating $v : $val"
      vnd.val .= val.val
      #calculate and fill in covariance
      #TODO rather broadcast than make new memory
      vnd.bw = val.cov
    end
  # elseif length(lsfPriors(csmc.cliqSubFg)) == 0 #FIXME
  #   @error "Par-3, clique $(csmc.cliq.index) failed to converge in upsolve, but ignoring since no priors" result
  else
    @error "Par-3, clique $(csmc.cliq.index) failed to converge in upsolve" result

    # propagate error to cleanly exit all cliques?
    beliefMsg = LikelihoodMessage(status=:ERROR_STATUS)
    for e in getEdgesParent(csmc.tree, csmc.cliq)
      @info "Par-3, $(csmc.cliq.index): generate up error on edge $(isa(e,Graphs.Edge) ? e.index : e)"
      putBeliefMessageUp!(csmc.tree, e, beliefMsg)
    end

    if length(getParent(csmc.tree, csmc.cliq)) == 0
      @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
        @info "Par-3 Root $(csmc.cliq.index): generate down error on edge $(isa(e,Graphs.Edge) ? e.index : e)"
        @async putBeliefMessageDown!(csmc.tree, e, LikelihoodMessage(status=:ERROR_STATUS))
      end
      @error "Par-3 $(csmc.cliq.index): Exit with error state"
      return IncrementalInference.exitStateMachine
    end

    return waitForDown_ParametricStateMachine
  end

  # Done with solve delete factors
  #TODO confirm, maybe don't delete mesage factors on subgraph, maybe delete if its priors, but not conditionals
  deleteMsgFactors!(csmc.cliqSubFg)

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_afterupsolve")

  #fill in belief
  #TODO createBeliefMessageParametric(csmc.cliqSubFg, csmc.cliq, solvekey=opts.solvekey)
  cliqSeparatorVarIds = getCliqSeparatorVarIds(csmc.cliq)
  #Fil in CliqueLikelihood
  cliqlikelihood = calculateMarginalCliqueLikelihood(vardict, Σ, varIds, cliqSeparatorVarIds)
  # @info "$(csmc.cliq.index) clique likelihood message $(cliqlikelihood)"
  beliefMsg = LikelihoodMessage(status=:UPSOLVED, variableOrder=cliqSeparatorVarIds, cliqueLikelihood=cliqlikelihood, msgType=ParametricMessage())

  #FIXME bit of a hack, only fill in variable beliefs if there are priors or for now more than one seperator
  if length(lsfPriors(csmc.cliqSubFg)) > 0 || length(cliqSeparatorVarIds) > 1
    for si in cliqSeparatorVarIds
      vnd = getSolverData(getVariable(csmc.cliqSubFg, si), :parametric)
      beliefMsg.belief[si] = TreeBelief(deepcopy(vnd))
    end
  end

  for e in getEdgesParent(csmc.tree, csmc.cliq)
    @info "$(csmc.cliq.index): put! on edge $(isa(e,Graphs.Edge) ? e.index : e)"
    putBeliefMessageUp!(csmc.tree, e, beliefMsg)
  end

  return waitForDown_ParametricStateMachine
end


"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 4
"""
function waitForDown_ParametricStateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "Par-4, wait for down messages if needed")
  setCliqDrawColor(csmc.cliq, "turquoise")

  for e in getEdgesParent(csmc.tree, csmc.cliq)
    @info "$(csmc.cliq.index): take! on edge $(isa(e,Graphs.Edge) ? e.index : e)"
    # Blocks until data is available.
    beliefMsg = takeBeliefMessageDown!(csmc.tree, e) # take!(csmc.tree.messages[e.index].downMsg)
    @info "$(csmc.cliq.index): Belief message received with status $(beliefMsg.status)"

    # logCSM(csmc, "Par-4 down msg"; beliefMsg)
    @debug "Par-4 down msg" beliefMsg
    # save down incoming message for use and debugging
    messages(csmc.cliq).downRx = beliefMsg
  
    # Down branching on message happens here
    
    # ERROR_STATUS
    if beliefMsg.status == :ERROR_STATUS

      setCliqDrawColor(csmc.cliq, "red")
      # csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

      @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
        @info "Par-4, $(csmc.cliq.index): put! error on edge $(isa(e,Graphs.Edge) ? e.index : e)"
        @async putBeliefMessageDown!(csmc.tree, e, LikelihoodMessage(status=:ERROR_STATUS))
      end
      @error "Par-4, $(csmc.cliq.index): Exit with error state"
      return IncrementalInference.exitStateMachine
    
    elseif beliefMsg.status == :DOWNSOLVED 
      return solveDown_ParametricStateMachine
    else
      error("waitForDown State Error: Unknown/unimplemented transision.")
    end
  end

  # The clique is a root
  # root clique down branching happens here
  return solveDown_ParametricStateMachine

end


"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 5
"""
function solveDown_ParametricStateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "Par-5, Solving down")

  setCliqDrawColor(csmc.cliq, "red")
  # csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  # TODO create function: 
  # updateMsgSeparators!(csmc.cliqSubFg, downmsg)
  downmsg = messages(csmc.cliq).downRx  #see #855
  svars = getCliqSeparatorVarIds(csmc.cliq)
  if !isnothing(downmsg)
    for (msym, belief) in downmsg.belief
      if msym in svars
        #TODO maybe combine variable and factor in new prior?
        vnd = getSolverData(getVariable(csmc.cliqSubFg, msym), :parametric)
        @info "$(csmc.cliq.index): Updating separator $msym from message $(belief.val)"
        vnd.val .= belief.val
        vnd.bw .= belief.bw
      end
    end
  end

  # store the cliqSubFg for later debugging
  # NOTE ITS not changed for now but keep here for possible future use
  # _dbgCSMSaveSubFG(csmc, "fg_beforedownsolve")

  # DownSolve cliqSubFg
  #only down solve if its not a root
  if length(getParent(csmc.tree, csmc.cliq)) != 0
    frontals = getCliqFrontalVarIds(csmc.cliq)
    vardict, result, flatvars, Σ = solveConditionalsParametric(csmc.cliqSubFg, frontals)
    #TEMP testing difference
    # vardict, result = solveFactorGraphParametric(csmc.cliqSubFg)
    # Pack all results in variables
    if result.g_converged || result.f_converged
      @info "$(csmc.cliq.index): subfg optim converged updating variables"
      for (v,val) in vardict
        @info "$(csmc.cliq.index) down: updating $v : $val"
        vnd = getSolverData(getVariable(csmc.cliqSubFg, v), :parametric)
        #Update subfg variables
        vnd.val .= val.val
        vnd.bw .= val.cov
      end
    else
      @error "Par-5, clique $(csmc.cliq.index) failed to converge in down solve" result

      #propagate error to cleanly exit all cliques?
      beliefMsg = LikelihoodMessage(status=:ERROR_STATUS)
      @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
        @info "Par-5, $(csmc.cliq.index): put! error on edge $(isa(e,Graphs.Edge) ? e.index : e)"
        @async putBeliefMessageDown!(csmc.tree, e, beliefMsg)#put!(csmc.tree.messages[e.index].downMsg, beliefMsg)
      end
      @error "Par-5, $(csmc.cliq.index): Exit with error state"
      return IncrementalInference.exitStateMachine

    end
  end

  #TODO fill in belief
  cliqFrontalVarIds = getCliqFrontalVarIds(csmc.cliq)
  #TODO createBeliefMessageParametric
  # beliefMsg = createBeliefMessageParametric(csmc.cliqSubFg, cliqFrontalVarIds, solvekey=opts.solvekey)
  beliefMsg = LikelihoodMessage(status=:DOWNSOLVED, msgType=ParametricMessage())
  for fi in cliqFrontalVarIds
    vnd = getSolverData(getVariable(csmc.cliqSubFg, fi), :parametric)
    beliefMsg.belief[fi] = TreeBelief(vnd)
    # beliefMsg.belief[fi] = TreeBelief(vnd.val, vnd.bw, vnd.inferdim, vnd.softtype.manifolds)
    @info "$(csmc.cliq.index): down message $fi : $beliefMsg"
  end

  # pass through the frontal variables that were sent from above
  if !isnothing(downmsg)
    pass_through_separators = intersect(svars, keys(downmsg.belief))
    for si in pass_through_separators
      beliefMsg.belief[si] = downmsg.belief[si]
      logCSM(csmc, "adding parent message"; sym=si, msg=downmsg.belief[si])
    end
  end

  #TODO sendBeliefMessageParametric(csmc, beliefMsg)
  #TODO maybe send a specific message to only the child that needs it
  @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
    @info "$(csmc.cliq.index): put! on edge $(isa(e,Graphs.Edge) ? e.index : e)"
    @async putBeliefMessageDown!(csmc.tree, e, beliefMsg)#put!(csmc.tree.messages[e.index].downMsg, beliefMsg)
  end

  @info "$(csmc.cliq.index): Solve completed"

  if isa(csmc.dfg, DFG.InMemoryDFGTypes)
    #TODO update frontal variables here directly
    frontsyms = getFrontals(csmc.cliq)
    transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frontsyms, updatePPE=false, solveKey=:parametric)

    #solve finished change color
    setCliqDrawColor(csmc.cliq, "lightblue")
    # csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

    @info "$(csmc.cliq.index): Finished"
    return IncrementalInference.exitStateMachine
  else
    #seems like a nice place to update remote variables here
    return updateRemote_ParametricStateMachine
  end
end

"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 6
"""
function updateRemote_ParametricStateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "Par-6, Updating Remote")
  #TODO update frontal variables remotely here
  #NOTE a new state is made to allow for coms retries and error traping kind of behaviour

  @info "$(csmc.cliq.index): Finish en klaar"
  return IncrementalInference.exitStateMachine

end
