
"""
    $SIGNATURES

EXPERIMENTAL: Init and start state machine.
"""
function initStartCliqStateMachine_X!(dfg::AbstractDFG,
                                       tree::AbstractBayesTree,
                                       cliq::TreeClique,
                                       cliqKey::Int;
                                       oldcliqdata::BayesTreeNodeData=BayesTreeNodeData(),
                                       verbose::Bool=false,
                                       drawtree::Bool=false,
                                       show::Bool=false,
                                       incremental::Bool=true,
                                       limititers::Int=20,
                                       upsolve::Bool=true,
                                       downsolve::Bool=true,
                                       recordhistory::Bool=false,
                                       delay::Bool=false,
                                       logger::SimpleLogger=SimpleLogger(Base.stdout),
                                       solve_progressbar=nothing,
                                       algorithm::Symbol=:default)

  # NOTE use tree and messages for operations involving children and parents
  # TODO deprecate children and prnt clique copies
  children = TreeClique[]
  prnt = TreeClique[]

  destType = dfg isa InMemoryDFGTypes ? typeof(dfg) : InMemDFGType

  csmc = CliqStateMachineContainer(dfg, initfg(destType, solverParams=getSolverParams(dfg)),
                                   tree, cliq,
                                   prnt, children,
                                   incremental, drawtree, downsolve, delay,
                                   getSolverParams(dfg), Dict{Symbol,String}(), oldcliqdata, logger, 
                                   cliqKey, algorithm) 

  nxt = buildCliqSubgraph_X_StateMachine

  statemachine = StateMachine{CliqStateMachineContainer}(next=nxt, name="cliq$(cliq.index)")


  # store statemachine and csmc in task
  if dfg.solverParams.dbg || recordhistory
    task_local_storage(:statemachine, statemachine)
    task_local_storage(:csmc, csmc)
  end

  logCSM(csmc, "Clique $(csmc.cliq.index) starting", loglevel=Logging.Info)

  while statemachine(csmc, verbose=verbose, iterlimit=limititers, recordhistory=recordhistory)
    !isnothing(solve_progressbar) && next!(solve_progressbar)
  end

  return statemachine.history

end

"""
    $SIGNATURES

Build a sub factor graph for clique variables from the larger factor graph.

Notes
- Parametric State machine function nr.1
"""
function buildCliqSubgraph_X_StateMachine(csmc::CliqStateMachineContainer)
  # build a local subgraph for inference operations
  syms = getCliqAllVarIds(csmc.cliq)

  logCSM(csmc, "X-1, build subgraph syms=$(syms)")

  frontsyms = getCliqFrontalVarIds(csmc.cliq)
  sepsyms = getCliqSeparatorVarIds(csmc.cliq)
  buildCliqSubgraph!(csmc.cliqSubFg, csmc.dfg, frontsyms, sepsyms)

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_build")

  # go to 2 wait for up
  return waitForUp_X_StateMachine
end

"""
    $SIGNATURES

Branching up state
Notes
- State machine function nr. 2
"""
function waitForUp_X_StateMachine(csmc::CliqStateMachineContainer)

  logCSM(csmc, "X-2, wait for up messages if needed")

  # setCliqDrawColor(csmc.cliq, "olive") #TODO don't know if this is correct color

  # JT empty upRx buffer to save messages, TODO It may be ok not to empty 
  beliefMessages = empty!(messages(csmc.cliq).upRx)

  # take! messages from edges
  @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
    @async begin
      thisEdge = isa(e,Graphs.Edge) ? e.index : e
      logCSM(csmc, "$(csmc.cliq.index): take! on edge $thisEdge")
      # Blocks until data is available. -- take! model
      beliefMsg = takeBeliefMessageUp!(csmc.tree, e)
      beliefMessages[thisEdge] = beliefMsg
      logCSM(csmc, "$(csmc.cliq.index): Belief message received with status $(beliefMsg.status)"; msgvars = keys(beliefMsg.belief))
    end
  end

  # get all statuses from messages
  all_child_status = map(msg -> msg.status, values(beliefMessages))
  
  # Main Branching happens here - all up messages received

  # If one up error is received propagate ERROR_STATUS 
  if :ERROR_STATUS in all_child_status

    putErrorUp(csmc)
    #if its a root, propagate error down
    #FIXME rather check if no parents with function (hasParents or isRoot)
    if length(getParent(csmc.tree, csmc.cliq)) == 0
      putErrorDown(csmc)
      return IncrementalInference.exitStateMachine
    end
    
    return waitForDown_X_StateMachine

  elseif csmc.algorithm == :parametric 
    !all(all_child_status .== :UPSOLVED) && error("#FIXME")
    return solveUp_ParametricStateMachine

  elseif true #TODO Currently all up goes through solveUp_X 
    return solveUp_X_StateMachine

  else
    error("waitForUp State Error: Unknown transision.")
  end
  
end



"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 3
"""
function solveUp_X_StateMachine(csmc::CliqStateMachineContainer)

  logCSM(csmc, "X-3, Solving Up")

  setCliqDrawColor(csmc.cliq, "red")
  opts = getSolverParams(csmc.dfg)

  # add message factors from upRx: cached messages taken from children saved in this clique
  addMsgFactors!(csmc.cliqSubFg, messages(csmc.cliq).upRx, UpwardPass)
  logCSM(csmc, "messages for up"; upmsg=lsf(csmc.cliqSubFg, tags=[:LIKELIHOODMESSAGE]))

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_beforeupsolve")

  #UPSOLVE here
  all_child_status = map(msg -> msg.status, values(messages(csmc.cliq).upRx))

  if all(all_child_status .== :UPSOLVED) && !areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq) 
    logCSM(csmc, "All children upsolved, not init, try init then upsolve"; c=csmc.cliqKey)
    varorder = getCliqVarInitOrderUp(csmc.cliqSubFg)
    someInit = cycleInitByVarOrder!(csmc.cliqSubFg, varorder, logger=csmc.logger)
  end

  if all(all_child_status .== :UPSOLVED) && areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq) 
    logCSM(csmc, "X-3 doing upSolve -- all initialized")

    __doCliqUpSolveInitialized!(csmc)
    
    solveStatus = :UPSOLVED
    
    converged_and_happy = true

  elseif !areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq)
    
        # FIXME experimental init to whatever is in frontals
        # should work if linear manifold
        # hardcoded off 
        linear_on_manifold = false
        init_for_differential = begin
          allvars = getVariables(csmc.cliqSubFg)
          any_init = any(isInitialized.(allvars))
          is_root = isempty(getEdgesParent(csmc.tree, csmc.cliq)) 
          logCSM(csmc, "init_for_differential: "; c=csmc.cliqKey, is_root=is_root, any_init=any_init)
          linear_on_manifold && !is_root && !any_init
        end
        
        if init_for_differential
          frontal_vars = getVariable.(csmc.cliqSubFg,  getCliqFrontalVarIds(csmc.cliq))
          filter!(!isInitialized, frontal_vars)
          foreach(fvar->getSolverData(fvar).initialized = true, frontal_vars)
          logCSM(csmc, "init_for_differential: "; c=csmc.cliqKey,lbl=getLabel.(frontal_vars))
        end
        ## END experimental

    setCliqDrawColor(csmc.cliq, "green")

    logCSM(csmc, "X-3, Trying up init -- all not initialized"; c=csmc.cliqKey)
     
    # structure for all up message densities computed during this initialization procedure.
    varorder = getCliqVarInitOrderUp(csmc.cliqSubFg)
    someInit = cycleInitByVarOrder!(csmc.cliqSubFg, varorder, logger=csmc.logger)
    # is clique fully upsolved or only partially?
    # print out the partial init status of all vars in clique
    printCliqInitPartialInfo(csmc.cliqSubFg, csmc.cliq, csmc.logger)
    logCSM(csmc, "X-3, solveUp_X try init -- someInit=$someInit, varorder=$varorder"; c=csmc.cliqKey)
  
    someInit ? setCliqDrawColor(csmc.cliq, "darkgreen") :  setCliqDrawColor(csmc.cliq, "lightgreen")

    solveStatus = someInit ? :INIT : :NO_INIT
    converged_and_happy = true

        ## FIXME init to whatever is in frontals
        # set frontals init back to false
        if init_for_differential #experimental_sommer_init_to_whatever_is_in_frontals
          foreach(fvar->getSolverData(fvar).initialized = false, frontal_vars)
          if someInit 
            solveStatus = :UPSOLVED
          end
        end
        ## END EXPERIMENTAL

  else
    setCliqDrawColor(csmc.cliq, "brown")
    logCSM(csmc, "X-3, we are initialized but children need to init, don't do anything")
    solveStatus = :INIT
    converged_and_happy = true
  end
  

  if converged_and_happy

  else # something went wrong propagate error
    @error "X-3, something wrong with solve up" 
    # propagate error to cleanly exit all cliques
    putErrorUp(csmc)
    if length(getParent(csmc.tree, csmc.cliq)) == 0
      putErrorDown(csmc)
      return IncrementalInference.exitStateMachine
    end

    return waitForDown_X_StateMachine
  end

  #fill in belief
  beliefMsg = prepCliqueMsgUpConsolidated(csmc.cliqSubFg, csmc.cliq, solveStatus, logger=csmc.logger)

  logCSM(csmc, "X-3, prepCliqueMsgUpConsolidated", msgon=keys(beliefMsg.belief), beliefMsg=beliefMsg)
  
  setCliqueStatus!(csmc.cliq, solveStatus)

  # Done with solve delete factors
  # remove msg factors that were added to the subfg
  tags_ = getSolverParams(csmc.cliqSubFg).useMsgLikelihoods ? [:UPWARD_COMMON;] : [:LIKELIHOODMESSAGE;]
  msgfcts= deleteMsgFactors!(csmc.cliqSubFg, tags_)
  logCSM(csmc, "8g, doCliqUpsSolveInit.! -- status = $(solveStatus), removing $(tags_) factors, length=$(length(msgfcts))")

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_afterupsolve")

  #propagate belief
  for e in getEdgesParent(csmc.tree, csmc.cliq)
    logCSM(csmc, "$(csmc.cliq.index): put! on edge $(isa(e,Graphs.Edge) ? e.index : e)")
    messages(csmc.cliq).upTx = deepcopy(beliefMsg)
    putBeliefMessageUp!(csmc.tree, e, beliefMsg)
  end

  return waitForDown_X_StateMachine
end


"""
    $SIGNATURES

Notes
- State machine function waitForDown nr. 4
"""
function waitForDown_X_StateMachine(csmc::CliqStateMachineContainer)

  logCSM(csmc, "X-4, wait for down messages if needed")

  # setCliqDrawColor(csmc.cliq, "lime")
 
  for e in getEdgesParent(csmc.tree, csmc.cliq)
    logCSM(csmc, "$(csmc.cliq.index): take! on edge $(isa(e,Graphs.Edge) ? e.index : e)")
    # Blocks until data is available.
    beliefMsg = takeBeliefMessageDown!(csmc.tree, e) # take!(csmc.tree.messages[e.index].downMsg)
    logCSM(csmc, "$(csmc.cliq.index): Belief message received with status $(beliefMsg.status)")

    logCSM(csmc, "X-4 down msg on $(keys(beliefMsg.belief))"; beliefMsg=beliefMsg)
    # save down incoming message for use and debugging
    messages(csmc.cliq).downRx = beliefMsg

    # Down branching happens here
    
    # ERROR_STATUS
    if beliefMsg.status == :ERROR_STATUS
      putErrorDown(csmc)
      return IncrementalInference.exitStateMachine

    elseif csmc.algorithm == :parametric
      beliefMsg.status != :DOWNSOLVED && error("#FIXME")
      return solveDown_ParametricStateMachine
    elseif beliefMsg.status == :DOWNSOLVED 
      return solveDown_X_StateMachine
    elseif beliefMsg.status == :INIT || beliefMsg.status == :NO_INIT
      return tryDownInit_X_StateMachine
    else
      logCSM(csmc, "Unknown state"; status=beliefMsg.status, loglevel=Logging.Error, c=csmc.cliqKey)
      error("waitForDown State Error: Unknown/unimplemented transision.")
    end
  end

  if csmc.algorithm == :parametric
    # The clique is a root
    # root clique down branching happens here
    return solveDown_ParametricStateMachine
  end

  # Special root case 
  #TODO improve
  solveStatus = getCliqueStatus(csmc.cliq)
  logCSM(csmc, "root case"; status=solveStatus, c=csmc.cliqKey)
  if solveStatus in [:INIT, :NO_INIT]
    return tryDownInit_X_StateMachine
  elseif solveStatus == :UPSOLVED
    setCliqueStatus!(csmc.cliq, :DOWNSOLVED)
    return solveDown_X_StateMachine
  else
    error("unknown status root $solveStatus")
  end


end

function tryDownInit_X_StateMachine(csmc::CliqStateMachineContainer)

  setCliqDrawColor(csmc.cliq, "olive")

  if length(getParent(csmc.tree, csmc.cliq)) != 0
    
    dwnmsgs = messages(csmc.cliq).downRx
    
    msgfcts = addMsgFactors!(csmc.cliqSubFg, dwnmsgs, DownwardPass)
    
    logCSM(csmc, "X-3, Trying Down init -- all not initialized")
        
    
    # structure for all up message densities computed during this initialization procedure.
    # XXX
    dwnkeys_ = lsf(csmc.cliqSubFg, tags=[:DOWNWARD_COMMON;]) .|> x->ls(csmc.cliqSubFg, x)[1]
    initorder = getCliqInitVarOrderDown(csmc.cliqSubFg, csmc.cliq, dwnkeys_)
    # initorder = getCliqVarInitOrderUp(csmc.tree, csmc.cliq)

    someInit = cycleInitByVarOrder!(csmc.cliqSubFg, initorder, logger=csmc.logger)
    # is clique fully upsolved or only partially?
    # print out the partial init status of all vars in clique
    printCliqInitPartialInfo(csmc.cliqSubFg, csmc.cliq, csmc.logger)
    logCSM(csmc, "8m, tryInitCliq_StateMachine -- someInit=$someInit, varorder=$initorder")

    solveStatus = someInit ? :INIT : :NO_INIT
    
    deleteMsgFactors!(csmc.cliqSubFg, msgfcts) # msgfcts # TODO, use tags=[:LIKELIHOODMESSAGE], see #760
    logCSM(csmc, "tryDownInit_X_StateMachine - removing factors, length=$(length(msgfcts))")

    someInit ? setCliqDrawColor(csmc.cliq, "seagreen") :  setCliqDrawColor(csmc.cliq, "khaki")

  else
    solveStatus = getCliqueStatus(csmc.cliq)
  end
  
  #fill in belief
  beliefMsg = CliqDownMessage(csmc, solveStatus)

  logCSM(csmc, "msg to send down"; beliefMsg=beliefMsg)
  # pass through the frontal variables that were sent from above
  downmsg = messages(csmc.cliq).downRx
  svars = getCliqSeparatorVarIds(csmc.cliq)
  if !isnothing(downmsg)
    pass_through_separators = intersect(svars, keys(downmsg.belief))
    for si in pass_through_separators
      beliefMsg.belief[si] = downmsg.belief[si]
      logCSM(csmc, "adding parent message"; sym=si, msg=downmsg.belief[si])
    end
  end

  #TODO maybe send a specific message to only the child that needs it
  @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
    logCSM(csmc, "$(csmc.cliq.index): put! on edge $(isa(e,Graphs.Edge) ? e.index : e)")
    @async putBeliefMessageDown!(csmc.tree, e, beliefMsg)#put!(csmc.tree.messages[e.index].downMsg, beliefMsg)
  end
  
  # detete all message factors to start clean
  deleteMsgFactors!(csmc.cliqSubFg) 


  return waitForUp_X_StateMachine
  
end



function CliqDownMessage(csmc::CliqStateMachineContainer, status=:DOWNSOLVED)

  #JT TODO maybe use Tx buffer
  newDwnMsgs = LikelihoodMessage(status=status)

  # create all messages from subfg
  for mk in getCliqFrontalVarIds(csmc.cliq)
    v = getVariable(csmc.cliqSubFg, mk)
    if isInitialized(v)
      newDwnMsgs.belief[mk] = TreeBelief(v)
    end
  end

  logCSM(csmc, "cliq $(csmc.cliq.index), CliqDownMessage, allkeys=$(keys(newDwnMsgs.belief))")
 
  return newDwnMsgs
end

"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 5
"""
function solveDown_X_StateMachine(csmc::CliqStateMachineContainer)

  logCSM(csmc, "X-5, Solving down")

  setCliqDrawColor(csmc.cliq, "maroon")

  # DownSolve cliqSubFg
  #only down solve if its not a root
  if length(getParent(csmc.tree, csmc.cliq)) != 0
    
    # TODO we can monitor the solve here to give it a timeout
    # add messages, do downsolve, remove messages
    logCSM(csmc, "11, doCliqDownSolve_StateMachine")
    
    #XXX
    # get down msg from Rx buffer (saved in take!)
    dwnmsgs = messages(csmc.cliq).downRx
    logCSM(csmc, "11, doCliqDownSolve_StateMachine -- dwnmsgs=$(collect(keys(dwnmsgs.belief)))")
  
    __doCliqDownSolve!(csmc, dwnmsgs)
    
    logCSM(csmc, "11, doCliqDownSolve_StateMachine -- finished with downGibbsCliqueDensity, now update csmc")

    # update clique subgraph with new status
    # setCliqDrawColor(csmc.cliq, "lightblue")

    # remove msg factors that were added to the subfg
    rmFcts = deleteMsgFactors!(csmc.cliqSubFg)
    logCSM(csmc, "11, doCliqDownSolve_StateMachine -- removing up message factors, length=$(length(rmFcts))")

    # store the cliqSubFg for later debugging
    _dbgCSMSaveSubFG(csmc, "fg_afterdownsolve")
    #XXX

    converged_and_happy = true

    if converged_and_happy

    else
      @error "X-5, clique $(csmc.cliq.index) failed in down solve"
      #propagate error to cleanly exit all cliques
      putErrorDown(csmc)
      return IncrementalInference.exitStateMachine
    end
  end

  #fill in belief
  beliefMsg = CliqDownMessage(csmc)

  if length(keys(beliefMsg.belief)) == 0
    logCSM(csmc, "Empty message on clique frontals"; loglevel=Logging.Error)
  end

  logCSM(csmc, "msg to send down on $(keys(beliefMsg.belief))"; beliefMsg=beliefMsg)
  # pass through the frontal variables that were sent from above
  downmsg = messages(csmc.cliq).downRx
  svars = getCliqSeparatorVarIds(csmc.cliq)
  if !isnothing(downmsg)
    pass_through_separators = intersect(svars, keys(downmsg.belief))
    for si in pass_through_separators
      beliefMsg.belief[si] = downmsg.belief[si]
      logCSM(csmc, "adding parent message"; sym=si, msg=downmsg.belief[si])
    end
  end

  #TODO maybe send a specific message to only the child that needs it
  @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
    logCSM(csmc, "$(csmc.cliq.index): put! on edge $(isa(e,Graphs.Edge) ? e.index : e)")
    @async putBeliefMessageDown!(csmc.tree, e, beliefMsg)#put!(csmc.tree.messages[e.index].downMsg, beliefMsg)
  end

  logCSM(csmc, "$(csmc.cliq.index): clique solve completed")

  if isa(csmc.dfg, DFG.InMemoryDFGTypes)
    #Update frontal variables here 

    # set PPE and solved for all frontals
    for sym in getCliqFrontalVarIds(csmc.cliq)
      # set PPE in cliqSubFg
      setVariablePosteriorEstimates!(csmc.cliqSubFg, sym)
      # set solved flag
      vari = getVariable(csmc.cliqSubFg, sym)
      setSolvedCount!(vari, getSolvedCount(vari, :default)+1, :default )
    end

    # transfer results to main factor graph
    frsyms = getCliqFrontalVarIds(csmc.cliq)
    logCSM(csmc, "11, finishingCliq -- going for transferUpdateSubGraph! on $frsyms")
    transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frsyms, csmc.logger, updatePPE=true)

    #solve finished change color
    setCliqDrawColor(csmc.cliq, "lightblue")
    # csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

    logCSM(csmc, "Clique $(csmc.cliq.index) finished", loglevel=Logging.Info)
    return IncrementalInference.exitStateMachine
  else
    #seems like a nice place to update remote variables here
    return updateRemote_ExpStateMachine
  end
end
