## =========================================================================================
## Initialization Functions -- 0
## =========================================================================================

"""
    $SIGNATURES

Init and start state machine.
"""
function initStartCliqStateMachine!(dfg::AbstractDFG,
                                    tree::AbstractBayesTree,
                                    cliq::TreeClique,
                                    timeout::Union{Nothing, <:Real}=nothing;
                                    oldcliqdata::BayesTreeNodeData=BayesTreeNodeData(),
                                    verbose::Bool=false,
                                    verbosefid=stdout,
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
                                   cliq.index, algorithm) 

  !upsolve && !downsolve && error("must attempt either up or down solve")
  # nxt = buildCliqSubgraph_StateMachine
  nxt = setCliqueRecycling_StateMachine

  csmiter_cb = getSolverParams(dfg).drawCSMIters ? ((st::StateMachine)->(cliq.attributes["xlabel"] = st.iter)) : ((st)->())

  statemachine = StateMachine{CliqStateMachineContainer}(next=nxt, name="cliq$(cliq.index)")


  # store statemachine and csmc in task
  if dfg.solverParams.dbg || recordhistory
    task_local_storage(:statemachine, statemachine)
    task_local_storage(:csmc, csmc)
  end

  logCSM(csmc, "Clique $(csmc.cliq.index) starting", loglevel=Logging.Debug)
  
  #TODO
  # timeout
  # verbosefid=verbosefid
  # injectDelayBefore=injectDelayBefore

  while statemachine(csmc, timeout; verbose=verbose, verbosefid=verbosefid, verboseXtra=getCliqueStatus(csmc.cliq), iterlimit=limititers, recordhistory=recordhistory, housekeeping_cb=csmiter_cb)
    !isnothing(solve_progressbar) && next!(solve_progressbar)
  end

  return CSMHistoryTuple.(statemachine.history)

end


"""
    $SIGNATURES

Recycle clique setup for later uses

Notes
- State machine function 0a
"""
function setCliqueRecycling_StateMachine(csmc::CliqStateMachineContainer)
  
  oldstatus = getCliqueStatus(csmc.oldcliqdata)
  
  # canCliqMargRecycle
  if areCliqVariablesAllMarginalized(csmc.dfg, csmc.cliq)
    getCliqueData(csmc.cliq).allmarginalized = true
    setCliqueStatus!(csmc.cliq, MARGINALIZED)

  # canCliqIncrRecycle
  # check if should be trying and can recycle clique computations
  elseif csmc.incremental && oldstatus == DOWNSOLVED  
    csmc.cliq.data.isCliqReused = true
    setCliqueStatus!(csmc.cliq, UPRECYCLED)
  end
  logCSM(csmc, "CSM-0a Recycling clique $(csmc.cliqKey) from $oldstatus"; 
              incremental=csmc.cliq.data.isCliqReused, 
              marginalized=getCliqueData(csmc.cliq).allmarginalized)

  return buildCliqSubgraph_StateMachine

end

"""
    $SIGNATURES

Build a sub factor graph for clique variables from the larger factor graph.

Notes
- State machine function 0b
"""
function buildCliqSubgraph_StateMachine(csmc::CliqStateMachineContainer)
  # build a local subgraph for inference operations
  syms = getCliqAllVarIds(csmc.cliq)

  logCSM(csmc, "CSM-0b build subgraph syms=$(syms)")

  frontsyms = getCliqFrontalVarIds(csmc.cliq)
  sepsyms = getCliqSeparatorVarIds(csmc.cliq)
  buildCliqSubgraph!(csmc.cliqSubFg, csmc.dfg, frontsyms, sepsyms)

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_build")

  # go to 2 wait for up
  return waitForUp_StateMachine
end

## =========================================================================================
## Wait for up -- 1
## =========================================================================================

"""
    $SIGNATURES

Branching up state
Notes
- State machine function 1
- Common state for handeling messages with take! approach
"""
function waitForUp_StateMachine(csmc::CliqStateMachineContainer)

  logCSM(csmc, "CSM-1 Wait for up messages if needed")

  # setCliqueDrawColor!(csmc.cliq, "olive") #TODO don't know if this is correct color

  # JT empty upRx buffer to save messages, TODO It may be ok not to empty 
  beliefMessages = empty!(getMessageBuffer(csmc.cliq).upRx)

  # take! messages from edges
  @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
    @async begin
      thisEdge = isa(e,Graphs.Edge) ? e.index : e
      logCSM(csmc, "CSM-1 $(csmc.cliq.index): take! on edge $thisEdge")
      # Blocks until data is available. -- take! model
      beliefMsg = takeBeliefMessageUp!(csmc.tree, e)
      beliefMessages[thisEdge] = beliefMsg
      logCSM(csmc, "CSM-1 $(csmc.cliq.index): Belief message received with status $(beliefMsg.status)"; msgvars = keys(beliefMsg.belief))
    end
  end

  # get all statuses from messages
  all_child_status = map(msg -> msg.status, values(beliefMessages))
  
  # Main Branching happens here - all up messages received

  # If one up error is received propagate ERROR_STATUS 
  if ERROR_STATUS in all_child_status

    putErrorUp(csmc)
    #if its a root, propagate error down
    #FIXME rather check if no parents with function (hasParents or isRoot)
    if length(getParent(csmc.tree, csmc.cliq)) == 0
      putErrorDown(csmc)
      return IncrementalInference.exitStateMachine
    end
    
    return waitForDown_StateMachine

  elseif csmc.algorithm == :parametric 
    !all(all_child_status .== UPSOLVED) && error("#FIXME")
    return solveUp_ParametricStateMachine

  elseif true #TODO Currently all up goes through solveUp 
    return preUpSolve_StateMachine

  else
    error("CSM-1 waitForUp State Error: Unknown transision.")
  end
  
end

## =========================================================================================
## Up functions -- 2
## =========================================================================================

"""
    $SIGNATURES

Notes
- State machine function 2a
"""
function preUpSolve_StateMachine(csmc::CliqStateMachineContainer)

  all_child_status = map(msg -> msg.status, values(getMessageBuffer(csmc.cliq).upRx))
  
  logCSM(csmc, "CSM-2a preUpSolve_StateMachine with child status"; all_child_status=all_child_status)

  #TODO perhaps don't add for MARGINALIZED 
  # always add messages in case its needed for downsolve (needed for differential)
  # add message factors from upRx: cached messages taken from children saved in this clique
  addMsgFactors!(csmc.cliqSubFg, getMessageBuffer(csmc.cliq).upRx, UpwardPass)
  logCSM(csmc, "CSM-2a messages for up"; upmsg=lsf(csmc.cliqSubFg, tags=[:LIKELIHOODMESSAGE]))

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_beforeupsolve")


  all_child_finished_up = all(in.(all_child_status, Ref([UPSOLVED, UPRECYCLED, MARGINALIZED])))

  #try to skip upsolve 
  if !getSolverParams(csmc.dfg).upsolve 
    return tryDownSolveOnly_StateMachine
  end

  #Clique and children UPSOLVED, UPRECYCLED or MARGINALIZED (finished upsolve)
  #no need to solve
  if getCliqueStatus(csmc.cliq) in [UPSOLVED, UPRECYCLED, MARGINALIZED] && all_child_finished_up
    logCSM(csmc, "CSM-2a Reusing clique $(csmc.cliqKey) as $(getCliqueStatus(csmc.cliq))")
    getCliqueStatus(csmc.cliq) == MARGINALIZED &&  setCliqueDrawColor!(csmc.cliq, "blue")
    return postUpSolve_StateMachine
  end

  # if all(all_child_status .== UPSOLVED) 
  if all_child_finished_up
    return solveUp_StateMachine
  elseif !areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq)
    return initUp_StateMachine
  else
    setCliqueDrawColor!(csmc.cliq, "brown")
    logCSM(csmc, "CSM-2a Clique $(csmc.cliqKey) is initialized but children need to init, don't do anything")
    setCliqueStatus!(csmc.cliq, INITIALIZED)
    return postUpSolve_StateMachine
  end
end

"""
  $SIGNATURES

Notes
- State machine function 2b
"""
function initUp_StateMachine(csmc)

        # FIXME experimental init to whatever is in frontals
        # should work if linear manifold
        # hardcoded off 
        linear_on_manifold = false
        init_for_differential = begin
          allvars = getVariables(csmc.cliqSubFg)
          any_init = any(isInitialized.(allvars))
          is_root = isempty(getEdgesParent(csmc.tree, csmc.cliq)) 
          logCSM(csmc, "CSM-2b init_for_differential: "; c=csmc.cliqKey, is_root=is_root, any_init=any_init)
          linear_on_manifold && !is_root && !any_init
        end
        
        if init_for_differential
          frontal_vars = getVariable.(csmc.cliqSubFg,  getCliqFrontalVarIds(csmc.cliq))
          filter!(!isInitialized, frontal_vars)
          foreach(fvar->getSolverData(fvar).initialized = true, frontal_vars)
          logCSM(csmc, "CSM-2b init_for_differential: "; c=csmc.cliqKey,lbl=getLabel.(frontal_vars))
        end
        ## END experimental

    setCliqueDrawColor!(csmc.cliq, "green")

    logCSM(csmc, "CSM-2b Trying up init -- all not initialized"; c=csmc.cliqKey)
     
    # structure for all up message densities computed during this initialization procedure.
    varorder = getCliqVarInitOrderUp(csmc.cliqSubFg)
    someInit = cycleInitByVarOrder!(csmc.cliqSubFg, varorder, logger=csmc.logger)
    # is clique fully upsolved or only partially?
    # print out the partial init status of all vars in clique
    printCliqInitPartialInfo(csmc.cliqSubFg, csmc.cliq, csmc.logger)
    logCSM(csmc, "CSM-2b solveUp try init -- someInit=$someInit, varorder=$varorder"; c=csmc.cliqKey)
  
    someInit ? setCliqueDrawColor!(csmc.cliq, "darkgreen") :  setCliqueDrawColor!(csmc.cliq, "lightgreen")

    solveStatus = someInit ? INITIALIZED : NO_INIT

        ## FIXME init to whatever is in frontals
        # set frontals init back to false
        if init_for_differential #experimental_sommer_init_to_whatever_is_in_frontals
          foreach(fvar->getSolverData(fvar).initialized = false, frontal_vars)
          if someInit 
            solveStatus = UPSOLVED
          end
        end
        ## END EXPERIMENTAL
  
  setCliqueStatus!(csmc.cliq, solveStatus)
  
  return postUpSolve_StateMachine

end


"""
  $SIGNATURES

Notes
- State machine function 2c
"""
function solveUp_StateMachine(csmc::CliqStateMachineContainer)
  
  logCSM(csmc, "CSM-2c Solving Up")

  setCliqueDrawColor!(csmc.cliq, "red")

  #Make sure all are initialized
  if !areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq) 
    logCSM(csmc, "CSM-2c All children upsolved, not init, try init then upsolve"; c=csmc.cliqKey)
    varorder = getCliqVarInitOrderUp(csmc.cliqSubFg)
    someInit = cycleInitByVarOrder!(csmc.cliqSubFg, varorder, logger=csmc.logger)
  end

  # Check again  
  if areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq) 
    logCSM(csmc, "CSM-2c doing upSolve -- all initialized")

    __doCliqUpSolveInitialized!(csmc)
    
    setCliqueStatus!(csmc.cliq, UPSOLVED)
  
  else
    _dbgCSMSaveSubFG(csmc, "fg_child_solved_cant_init")
    # it can be a leaf
    logCSM(csmc, "CSM-2c solveUp -- all children upsolved, but init failed.")
  end
  
  # if converged_and_happy

  # else # something went wrong propagate error
  #   @error "X-3, something wrong with solve up" 
  #   # propagate error to cleanly exit all cliques
  #   putErrorUp(csmc)
  #   if length(getParent(csmc.tree, csmc.cliq)) == 0
  #     putErrorDown(csmc)
  #     return IncrementalInference.exitStateMachine
  #   end

  #   return waitForDown_StateMachine
  # end


  return postUpSolve_StateMachine
end

"""
  $SIGNATURES

CSM function only called when `getSolverParams(dfg).upsolve == false` that tries to skip upsolve.
Notes
- Cliques are uprecycled to add differential messages. 
- State machine function 2d
"""
function tryDownSolveOnly_StateMachine(csmc::CliqStateMachineContainer)
  logCSM(csmc, "CSM-2d tryDownSolveOnly_StateMachine clique $(csmc.cliqKey) status $(getCliqueStatus(csmc.cliq))")

  logCSM(csmc, "CSM-2d Skipping upsolve clique $(csmc.cliqKey)"; loglevel=Logging.Info, st=getCliqueStatus(csmc.cliq))
  if getCliqueStatus(csmc.cliq) == NULL 
    logCSM(csmc, "CSM-2d Clique $(csmc.cliqKey) status NULL, trying as UPRECYCLED"; loglevel=Logging.Warn)
    
    # Are all variables solved at least once?
    if all(getSolvedCount.(getVariables(csmc.cliqSubFg)) .> 0)
      setCliqueStatus!(csmc.cliq, UPRECYCLED)
    else
      logCSM(csmc, "CSM-2d Clique $(csmc.cliqKey) cannot be UPRECYCLED, all variables not solved. Set solverParams to upsolve=true.";
             loglevel=Logging.Error)
      # propagate error to cleanly exit all cliques
      putErrorUp(csmc)
      if length(getParent(csmc.tree, csmc.cliq)) == 0
        putErrorDown(csmc)
        return IncrementalInference.exitStateMachine
      end
      return waitForDown_StateMachine

    end

  end

  return postUpSolve_StateMachine
end

"""
    $SIGNATURES

Post-upsolve remove message factors and send messages
Notes
- State machine function 2e
"""
function postUpSolve_StateMachine(csmc::CliqStateMachineContainer)

  solveStatus = getCliqueStatus(csmc.cliq)
  #fill in belief
  beliefMsg = prepCliqueMsgUpConsolidated(csmc.cliqSubFg, csmc.cliq, solveStatus, logger=csmc.logger)

  logCSM(csmc, "CSM-2e prepCliqueMsgUpConsolidated", msgon=keys(beliefMsg.belief), beliefMsg=beliefMsg)

  # Done with solve delete factors
  # remove msg factors that were added to the subfg
  tags_ = getSolverParams(csmc.cliqSubFg).useMsgLikelihoods ? [:UPWARD_COMMON;] : [:LIKELIHOODMESSAGE;]
  msgfcts= deleteMsgFactors!(csmc.cliqSubFg, tags_)
  logCSM(csmc, "CSM-2e doCliqUpsSolveInit.! -- status = $(solveStatus), removing $(tags_) factors, length=$(length(msgfcts))")

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_afterupsolve")

  #propagate belief
  for e in getEdgesParent(csmc.tree, csmc.cliq)
    logCSM(csmc, "CSM-2e $(csmc.cliq.index): put! on edge $(isa(e,Graphs.Edge) ? e.index : e)")
    getMessageBuffer(csmc.cliq).upTx = deepcopy(beliefMsg)
    putBeliefMessageUp!(csmc.tree, e, beliefMsg)
  end

  if getSolverParams(csmc.dfg).downsolve
    return waitForDown_StateMachine
  else
    return updateFromSubgraph_StateMachine
  end
end

## =========================================================================================
## Wait for Down -- 3
## =========================================================================================

"""
    $SIGNATURES

Notes
- State machine function waitForDown 3
"""
function waitForDown_StateMachine(csmc::CliqStateMachineContainer)

  logCSM(csmc, "CSM-3 wait for down messages if needed")

  # setCliqueDrawColor!(csmc.cliq, "lime")
 
  for e in getEdgesParent(csmc.tree, csmc.cliq)
    logCSM(csmc, "CSM-3 $(csmc.cliq.index): take! on edge $(isa(e,Graphs.Edge) ? e.index : e)")
    # Blocks until data is available.
    beliefMsg = takeBeliefMessageDown!(csmc.tree, e) # take!(csmc.tree.messageChannels[e.index].downMsg)
    logCSM(csmc, "CSM-3 $(csmc.cliq.index): Belief message received with status $(beliefMsg.status)")

    logCSM(csmc, "CSM-3 down msg on $(keys(beliefMsg.belief))"; beliefMsg=beliefMsg)
    # save down incoming message for use and debugging
    getMessageBuffer(csmc.cliq).downRx = beliefMsg

    # Down branching happens here
    
    # ERROR_STATUS
    if beliefMsg.status == ERROR_STATUS
      putErrorDown(csmc)
      return IncrementalInference.exitStateMachine

    elseif csmc.algorithm == :parametric
      beliefMsg.status != DOWNSOLVED && error("#FIXME")
      return solveDown_ParametricStateMachine
    elseif beliefMsg.status in [MARGINALIZED, DOWNSOLVED, INITIALIZED, NO_INIT]
      return preDownSolve_StateMachine
    # elseif beliefMsg.status == DOWNSOLVED 
    #   return solveDown_StateMachine
    # elseif beliefMsg.status == INITIALIZED || beliefMsg.status == NO_INIT
    #   return tryDownInit_StateMachine
    else
      logCSM(csmc, "CSM-3 Unknown state"; status=beliefMsg.status, loglevel=Logging.Error, c=csmc.cliqKey)
      error("CSM-3 waitForDown State Error: Unknown/unimplemented transision.")
    end
  end

  # The clique is a root
  # root clique down branching happens here
  if csmc.algorithm == :parametric
    return solveDown_ParametricStateMachine
  else 
    return preDownSolve_StateMachine
  end

end

## =========================================================================================
## Down Functions -- 4
## =========================================================================================

## TODO Consolidate
function CliqDownMessage(csmc::CliqStateMachineContainer, status=DOWNSOLVED)

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
- State machine function 4a
"""
function preDownSolve_StateMachine(csmc::CliqStateMachineContainer)
  logCSM(csmc, "CSM-4a Preparing for down init/solve")
  
  opts = getSolverParams(csmc.dfg)
  
  # get down msg from Rx buffer (saved in take!)
  dwnmsgs = getMessageBuffer(csmc.cliq).downRx
  
  # DownSolve cliqSubFg
  #only down solve if its not a root and not MARGINALIZED
  if length(getParent(csmc.tree, csmc.cliq)) != 0 && getCliqueStatus(csmc.cliq) != MARGINALIZED
    
    logCSM(csmc, "CSM-4a doCliqDownSolve_StateMachine -- dwnmsgs=$(collect(keys(dwnmsgs.belief)))")
    # maybe cycle through separators (or better yet, just use values directly -- see next line)
    msgfcts = addMsgFactors!(csmc.cliqSubFg, dwnmsgs, DownwardPass)
    
    if dwnmsgs.status in [DOWNSOLVED, MARGINALIZED] 
      logCSM(csmc, "CSM-4a doCliqDownSolve_StateMachine")
      return solveDown_StateMachine
    elseif dwnmsgs.status == INITIALIZED || dwnmsgs.status == NO_INIT
      return tryDownInit_StateMachine
    else
      logCSM(csmc, "CSM-4a Unknown state"; status=dwnmsgs.status, loglevel=Logging.Error, c=csmc.cliqKey)
      error("CSM-4a waitForDown State Error: Unknown/unimplemented transision.")
    end
  else 
    # Special root case or MARGINALIZED
    #TODO improve
    solveStatus = getCliqueStatus(csmc.cliq)
    logCSM(csmc, "CSM-4a root case or MARGINALIZED"; status=solveStatus, c=csmc.cliqKey)
    if solveStatus in [INITIALIZED, NO_INIT, UPSOLVED, UPRECYCLED, MARGINALIZED]
      solveStatus == MARGINALIZED &&  setCliqueDrawColor!(csmc.cliq, "blue")
      if solveStatus in [UPSOLVED, UPRECYCLED]
        setCliqueStatus!(csmc.cliq, DOWNSOLVED)
      end
      return postDownSolve_StateMachine
    else
      error("CSM-4a unknown status root $solveStatus")
    end
  end

  
end

"""
    $SIGNATURES

Notes
- State machine function 4b
"""
function tryDownInit_StateMachine(csmc::CliqStateMachineContainer)

  setCliqueDrawColor!(csmc.cliq, "olive")

  logCSM(csmc, "CSM-4b Trying Down init -- all not initialized") 
    
  # structure for all up message densities computed during this initialization procedure.
  # XXX
  dwnkeys_ = lsf(csmc.cliqSubFg, tags=[:DOWNWARD_COMMON;]) .|> x->ls(csmc.cliqSubFg, x)[1]
  initorder = getCliqInitVarOrderDown(csmc.cliqSubFg, csmc.cliq, dwnkeys_)
  # initorder = getCliqVarInitOrderUp(csmc.tree, csmc.cliq)

  someInit = cycleInitByVarOrder!(csmc.cliqSubFg, initorder, logger=csmc.logger)
  # is clique fully upsolved or only partially?
  # print out the partial init status of all vars in clique
  printCliqInitPartialInfo(csmc.cliqSubFg, csmc.cliq, csmc.logger)
  logCSM(csmc, "CSM-4b tryInitCliq_StateMachine -- someInit=$someInit, varorder=$initorder")

  
  msgfcts = deleteMsgFactors!(csmc.cliqSubFg, [:DOWNWARD_COMMON;]) # msgfcts # TODO, use tags=[:LIKELIHOODMESSAGE], see #760
  logCSM(csmc, "CSM-4b tryDownInit_StateMachine - removing factors, length=$(length(msgfcts))")
  
  solveStatus = someInit ? INITIALIZED : NO_INIT
  someInit ? setCliqueDrawColor!(csmc.cliq, "seagreen") :  setCliqueDrawColor!(csmc.cliq, "khaki")
  setCliqueStatus!(csmc.cliq, solveStatus)

  return postDownSolve_StateMachine
end


"""
    $SIGNATURES

Notes
- State machine function 4c
"""
function solveDown_StateMachine(csmc::CliqStateMachineContainer)

  logCSM(csmc, "CSM-4c Solving down")

  setCliqueDrawColor!(csmc.cliq, "maroon")

  # DownSolve cliqSubFg
  # add messages, do downsolve, remove messages
  
  #XXX
  # get down msg from Rx buffer (saved in take!)
  dwnmsgs = getMessageBuffer(csmc.cliq).downRx

  opts = getSolverParams(csmc.dfg)

  # maybe cycle through separators (or better yet, just use values directly -- see next line)
  msgfcts = addMsgFactors!(csmc.cliqSubFg, dwnmsgs, DownwardPass)

  # force separator variables in cliqSubFg to adopt down message values
  updateSubFgFromDownMsgs!(csmc.cliqSubFg, dwnmsgs, getCliqSeparatorVarIds(csmc.cliq))

  #XXX test with and without
  # add required all frontal connected factors
  if !opts.useMsgLikelihoods
    newvars, newfcts = addDownVariableFactors!(csmc.dfg, csmc.cliqSubFg, csmc.cliq, csmc.logger, solvable=1)
  end

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_beforedownsolve")
  

  ## new way
  # calculate belief on each of the frontal variables and iterate if required
  solveCliqDownFrontalProducts!(csmc.cliqSubFg, csmc.cliq, opts, csmc.logger)
  csmc.dodownsolve = false
  
  logCSM(csmc, "CSM-4c solveDown -- finished with downGibbsCliqueDensity, now update csmc")

  # update clique subgraph with new status
  # setCliqueDrawColor!(csmc.cliq, "lightblue")

  # remove msg factors that were added to the subfg
  rmFcts = deleteMsgFactors!(csmc.cliqSubFg)
  logCSM(csmc, "CSM-4c solveDown -- removing up message factors, length=$(length(rmFcts))")

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_afterdownsolve")

  setCliqueStatus!(csmc.cliq, DOWNSOLVED) 
  
  logCSM(csmc, "CSM-4c $(csmc.cliq.index): clique down solve completed")

  return postDownSolve_StateMachine

end


"""
    $SIGNATURES

Notes
- State machine function 4d
"""
function postDownSolve_StateMachine(csmc::CliqStateMachineContainer)
  
  solveStatus = getCliqueStatus(csmc.cliq)
  #fill in belief
  #TODO use prepSetCliqueMsgDownConsolidated
  beliefMsg = CliqDownMessage(csmc, solveStatus)

  if length(keys(beliefMsg.belief)) == 0
    logCSM(csmc, "CSM-4d Empty message on clique $(csmc.cliqKey) frontals"; loglevel=Logging.Info)
  end

  logCSM(csmc, "CSM-4d msg to send down on $(keys(beliefMsg.belief))"; beliefMsg=beliefMsg)
  # pass through the frontal variables that were sent from above
  downmsg = getMessageBuffer(csmc.cliq).downRx
  svars = getCliqSeparatorVarIds(csmc.cliq)
  if !isnothing(downmsg)
    pass_through_separators = intersect(svars, keys(downmsg.belief))
    for si in pass_through_separators
      beliefMsg.belief[si] = downmsg.belief[si]
      logCSM(csmc, "CSM-4d adding parent message"; sym=si, msg=downmsg.belief[si])
    end
  end

  #TODO maybe send a specific message to only the child that needs it
  @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
    logCSM(csmc, "CSM-4d $(csmc.cliq.index): put! on edge $(isa(e,Graphs.Edge) ? e.index : e)")
    @async putBeliefMessageDown!(csmc.tree, e, beliefMsg)#put!(csmc.tree.messageChannels[e.index].downMsg, beliefMsg)
  end
  
  if getCliqueStatus(csmc.cliq) in [DOWNSOLVED, MARGINALIZED]

    return updateFromSubgraph_StateMachine
  
  else
    # detete all message factors to start clean
    deleteMsgFactors!(csmc.cliqSubFg) 

    return waitForUp_StateMachine
  end

end

## =========================================================================================
## Finalize Functions -- 5
## =========================================================================================

"""
    $SIGNATURES

The last step in CSM to update the main FG from the sub FG.

Notes
- CSM function 5
"""
function updateFromSubgraph_StateMachine(csmc::CliqStateMachineContainer)
  
  # NOTE possible future use for things like retry on CGDFGs 
  # if isa(csmc.dfg, DFG.InMemoryDFGTypes)
  # else
  #   #seems like a nice place to update remote variables here
  #   return updateRemote_ExpStateMachine
  # end

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
  logCSM(csmc, "CSM-5 finishingCliq -- going for transferUpdateSubGraph! on $frsyms")
  transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frsyms, csmc.logger, updatePPE=true)

  #solve finished change color
  setCliqueDrawColor!(csmc.cliq, "lightblue")

  logCSM(csmc, "CSM-5 Clique $(csmc.cliq.index) finished", loglevel=Logging.Info)
  return IncrementalInference.exitStateMachine

end
