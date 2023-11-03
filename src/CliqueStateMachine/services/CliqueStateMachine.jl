## =========================================================================================
## Initialization Functions -- 0
## =========================================================================================

"""
    $SIGNATURES

Init and start state machine.
"""
function initStartCliqStateMachine!(
  dfg::AbstractDFG,
  tree::AbstractBayesTree,
  cliq::TreeClique,
  timeout::Union{Nothing, <:Real} = nothing;
  oldcliqdata::BayesTreeNodeData = BayesTreeNodeData(),
  verbose::Bool = false,
  verbosefid = stdout,
  drawtree::Bool = false,
  show::Bool = false,
  incremental::Bool = true,
  limititers::Int = 20,
  upsolve::Bool = true,
  downsolve::Bool = true,
  recordhistory::Bool = false,
  delay::Bool = false,
  logger::SimpleLogger = SimpleLogger(Base.stdout),
  solve_progressbar = nothing,
  algorithm::Symbol = :default,
  solveKey::Symbol = algorithm,
)

  # NOTE use tree and messages for operations involving children and parents
  # TODO deprecate children and prnt clique copies
  # children = TreeClique[]
  # prnt = TreeClique[]

  destType = dfg isa InMemoryDFGTypes ? typeof(dfg) : LocalDFG

  csmc = CliqStateMachineContainer(
    dfg,
    initfg(destType; solverParams = getSolverParams(dfg)),
    tree,
    cliq,
    incremental,
    drawtree,
    downsolve,
    delay,
    getSolverParams(dfg),
    Dict{Symbol, String}(),
    oldcliqdata,
    logger,
    cliq.id,
    algorithm,
    0,
    true,
    solveKey,
    0,
  )

  !upsolve && !downsolve && error("must attempt either up or down solve")
  # nxt = buildCliqSubgraph_StateMachine
  nxt = setCliqueRecycling_StateMachine

  csmiter_cb = if getSolverParams(dfg).drawCSMIters
    ((st::StateMachine) -> (cliq.attributes["xlabel"] = st.iter; csmc._csm_iter = st.iter))
  else
    ((st) -> (csmc._csm_iter = st.iter))
  end

  statemachine =
    StateMachine{CliqStateMachineContainer}(; next = nxt, name = "cliq$(getId(cliq))")

  # store statemachine and csmc in task
  if dfg.solverParams.dbg || recordhistory
    task_local_storage(:statemachine, statemachine)
    task_local_storage(:csmc, csmc)
  end

  logCSM(csmc, "Clique $(getId(csmc.cliq)) starting"; loglevel = Logging.Debug)

  #TODO
  # timeout
  # verbosefid=verbosefid
  # injectDelayBefore=injectDelayBefore

  while statemachine(
    csmc,
    timeout;
    verbose = verbose,
    verbosefid = verbosefid,
    verboseXtra = getCliqueStatus(csmc.cliq),
    iterlimit = limititers,
    recordhistory = recordhistory,
    housekeeping_cb = csmiter_cb,
  )
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
  logCSM(
    csmc,
    "CSM-0a Recycling clique $(csmc.cliqId) from $oldstatus";
    incremental = csmc.cliq.data.isCliqReused,
    marginalized = getCliqueData(csmc.cliq).allmarginalized,
  )

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

  # TODO optimize by only fetching csmc.solveKey -- upgrades required
  buildCliqSubgraph!(csmc.cliqSubFg, csmc.dfg, frontsyms, sepsyms)

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_build")

  # go to 2 wait for up
  return presolveChecklist_StateMachine
end

"""
    $SIGNATURES

Check that the csmc container has everything it needs to proceed with init-ference.

DevNotes
- TODO marginalized flag might be wrong default.
"""
function presolveChecklist_StateMachine(csmc::CliqStateMachineContainer)

  # check if solveKey is available in all variables?
  for var in getVariable.(csmc.cliqSubFg, ls(csmc.cliqSubFg))
    if !(csmc.solveKey in listSolveKeys(var))
      logCSM(
        csmc,
        "CSM-0b create empty data for $(getLabel(var)) on solveKey=$(csmc.solveKey)",
      )
      varType = getVariableType(var)
      # FIXME check the marginalization requirements
      setDefaultNodeData!(
        var,
        0,
        getSolverParams(csmc.cliqSubFg).N,
        getDimension(varType);
        solveKey = csmc.solveKey,
        initialized = false,
        varType = varType,
        dontmargin = false,
      )
      #
      @info "create vnd solveKey" csmc.solveKey N
      @info "also" listSolveKeys(var)
    end
  end

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
      thisEdge = e.dst
      logCSM(csmc, "CSM-1 $(csmc.cliq.id): take! on edge $thisEdge")
      # Blocks until data is available. -- take! model
      beliefMsg = takeBeliefMessageUp!(csmc.tree, e)
      beliefMessages[thisEdge] = beliefMsg
      logCSM(
        csmc,
        "CSM-1 $(csmc.cliq.id): Belief message received with status $(beliefMsg.status)";
        msgvars = keys(beliefMsg.belief),
      )
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

  logCSM(
    csmc,
    "CSM-2a preUpSolve_StateMachine with child status";
    all_child_status = all_child_status,
  )

  #TODO perhaps don't add for MARGINALIZED 
  # always add messages in case its needed for downsolve (needed for differential)
  # add message factors from upRx: cached messages taken from children saved in this clique
  addMsgFactors!(csmc.cliqSubFg, getMessageBuffer(csmc.cliq).upRx, UpwardPass)
  logCSM(
    csmc,
    "CSM-2a messages for up";
    upmsg = lsf(csmc.cliqSubFg; tags = [:__LIKELIHOODMESSAGE__]),
  )

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_beforeupsolve")

  all_child_finished_up =
    all(in.(all_child_status, Ref([UPSOLVED, UPRECYCLED, MARGINALIZED])))
  logCSM(
    csmc,
    "CSM-2a, clique $(csmc.cliqId) all_child_finished_up $(all_child_finished_up)",
  )

  #try to skip upsolve 
  if !getSolverParams(csmc.dfg).upsolve
    return tryDownSolveOnly_StateMachine
  end

  #Clique and children UPSOLVED, UPRECYCLED or MARGINALIZED (finished upsolve)
  #no need to solve
  if getCliqueStatus(csmc.cliq) in [UPSOLVED, UPRECYCLED, MARGINALIZED] &&
     all_child_finished_up
    logCSM(csmc, "CSM-2a Reusing clique $(csmc.cliqId) as $(getCliqueStatus(csmc.cliq))")
    getCliqueStatus(csmc.cliq) == MARGINALIZED && setCliqueDrawColor!(csmc.cliq, "blue")
    getCliqueStatus(csmc.cliq) == UPRECYCLED && setCliqueDrawColor!(csmc.cliq, "orange")
    return postUpSolve_StateMachine
  end

  # if all(all_child_status .== UPSOLVED) 
  if all_child_finished_up
    return solveUp_StateMachine
  elseif !areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq, csmc.solveKey)
    return initUp_StateMachine
  else
    setCliqueDrawColor!(csmc.cliq, "brown")
    logCSM(
      csmc,
      "CSM-2a Clique $(csmc.cliqId) is initialized but children need to init, don't do anything",
    )
    setCliqueStatus!(csmc.cliq, INITIALIZED)
    return postUpSolve_StateMachine
  end
end

"""
  $SIGNATURES

Notes
- State machine function 2b
"""
function initUp_StateMachine(csmc::CliqStateMachineContainer)
  csmc.init_iter += 1

  # FIXME experimental init to whatever is in frontals
  # should work if linear manifold
  # hardcoded off 
  linear_on_manifold = false
  init_for_differential = begin
    allvars = getVariables(csmc.cliqSubFg)
    any_init = any(isInitialized.(allvars, csmc.solveKey))
    is_root = isempty(getEdgesParent(csmc.tree, csmc.cliq))
    logCSM(
      csmc,
      "CSM-2b init_for_differential: ";
      c = csmc.cliqId,
      is_root = is_root,
      any_init = any_init,
    )
    linear_on_manifold && !is_root && !any_init
  end

  if init_for_differential
    frontal_vars = getVariable.(csmc.cliqSubFg, getCliqFrontalVarIds(csmc.cliq))
    filter!(!isInitialized, frontal_vars)
    foreach(fvar -> getSolverData(fvar, csmc.solveKey).initialized = true, frontal_vars)
    logCSM(
      csmc,
      "CSM-2b init_for_differential: ";
      c = csmc.cliqId,
      lbl = getLabel.(frontal_vars),
    )
  end
  ## END experimental

  setCliqueDrawColor!(csmc.cliq, "green")

  logCSM(csmc, "CSM-2b Trying up init -- all not initialized"; c = csmc.cliqId)

  # structure for all up message densities computed during this initialization procedure.
  varorder = getCliqVarInitOrderUp(csmc.cliqSubFg)
  someInit = cycleInitByVarOrder!(
    csmc.cliqSubFg,
    varorder;
    solveKey = csmc.solveKey,
    logger = csmc.logger,
  )
  # is clique fully upsolved or only partially?
  # print out the partial init status of all vars in clique
  printCliqInitPartialInfo(csmc.cliqSubFg, csmc.cliq, csmc.solveKey, csmc.logger)
  logCSM(
    csmc,
    "CSM-2b solveUp try init -- someInit=$someInit, varorder=$varorder";
    c = csmc.cliqId,
  )

  if someInit
    setCliqueDrawColor!(csmc.cliq, "darkgreen")
  else
    setCliqueDrawColor!(csmc.cliq, "lightgreen")
  end

  solveStatus = someInit ? INITIALIZED : NO_INIT

  ## FIXME init to whatever is in frontals
  # set frontals init back to false
  if init_for_differential #experimental_sommer_init_to_whatever_is_in_frontals
    foreach(fvar -> getSolverData(fvar, csmc.solveKey).initialized = false, frontal_vars)
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
  logCSM(csmc, "CSM-2c, cliq=$(csmc.cliqId) Solving Up")

  setCliqueDrawColor!(csmc.cliq, "red")

  #Make sure all are initialized
  if !areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq, csmc.solveKey)
    logCSM(
      csmc,
      "CSM-2c All children upsolved, not init, try init then upsolve";
      c = csmc.cliqId,
    )
    varorder = getCliqVarInitOrderUp(csmc.cliqSubFg)
    someInit = cycleInitByVarOrder!(
      csmc.cliqSubFg,
      varorder;
      solveKey = csmc.solveKey,
      logger = csmc.logger,
    )
  end

  isinit = areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq, csmc.solveKey)
  logCSM(csmc, "CSM-2c midway, isinit=$isinit")
  # Check again  
  if isinit
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
  logCSM(
    csmc,
    "CSM-2d tryDownSolveOnly_StateMachine clique $(csmc.cliqId) status $(getCliqueStatus(csmc.cliq))",
  )

  logCSM(
    csmc,
    "CSM-2d Skipping upsolve clique $(csmc.cliqId)";
    loglevel = Logging.Info,
    st = getCliqueStatus(csmc.cliq),
  )
  if getCliqueStatus(csmc.cliq) == NULL
    logCSM(
      csmc,
      "CSM-2d Clique $(csmc.cliqId) status NULL, trying as UPRECYCLED";
      loglevel = Logging.Warn,
    )

    # Are all variables solved at least once?
    if all(getSolvedCount.(getVariables(csmc.cliqSubFg)) .> 0)
      setCliqueStatus!(csmc.cliq, UPRECYCLED)
    else
      logCSM(
        csmc,
        "CSM-2d Clique $(csmc.cliqId) cannot be UPRECYCLED, all variables not solved. Set solverParams to upsolve=true.";
        loglevel = Logging.Error,
      )
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
  # fill in belief
  logCSM(csmc, "CSM-2e prepCliqueMsgUp, going for prepCliqueMsgUp")
  beliefMsg = prepCliqueMsgUp(
    csmc.cliqSubFg,
    csmc.cliq,
    csmc.solveKey,
    solveStatus;
    logger = csmc.logger,
    sender = (; id = csmc.cliq.id.value, step = csmc._csm_iter),
  )
  #

  logCSM(
    csmc,
    "CSM-2e prepCliqueMsgUp";
    msgon = keys(beliefMsg.belief),
    beliefMsg = beliefMsg,
  )

  # Done with solve delete factors
  # remove msg factors that were added to the subfg
  tags_ = if getSolverParams(csmc.cliqSubFg).useMsgLikelihoods
    [:__UPWARD_COMMON__;]
  else
    [:__LIKELIHOODMESSAGE__;]
  end
  msgfcts = deleteMsgFactors!(csmc.cliqSubFg, tags_)
  logCSM(
    csmc,
    "CSM-2e doCliqUpsSolveInit.! -- status = $(solveStatus), removing $(tags_) factors, length=$(length(msgfcts))",
  )

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_afterupsolve")

  # warn and clean exit on stalled tree init
  if csmc.init_iter > getSolverParams(csmc.cliqSubFg).limittreeinit_iters
    logCSM(
      csmc,
      "CSM-2e Clique $(csmc.cliqId) tree init failed, max init retries reached.";
      loglevel = Logging.Error,
    )
    putErrorUp(csmc)
    if length(getParent(csmc.tree, csmc.cliq)) == 0
      putErrorDown(csmc)
      return IncrementalInference.exitStateMachine
    end
    return waitForDown_StateMachine
  end

  # always put up belief message in upTx, only used for debugging isolated cliques
  getMessageBuffer(csmc.cliq).upTx = deepcopy(beliefMsg)
  #propagate belief
  for e in getEdgesParent(csmc.tree, csmc.cliq)
    logCSM(csmc, "CSM-2e $(csmc.cliq.id): put! on edge $(e)")
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
    logCSM(csmc, "CSM-3 $(csmc.cliq.id): take! on edge $(e)")
    # Blocks until data is available.
    beliefMsg = takeBeliefMessageDown!(csmc.tree, e) # take!(csmc.tree.messageChannels[e.index].downMsg)
    logCSM(
      csmc,
      "CSM-3 $(csmc.cliq.id): Belief message received with status $(beliefMsg.status)",
    )

    logCSM(csmc, "CSM-3 down msg on $(keys(beliefMsg.belief))"; beliefMsg = beliefMsg)
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
      logCSM(
        csmc,
        "CSM-3 Unknown state";
        status = beliefMsg.status,
        loglevel = Logging.Error,
        c = csmc.cliqId,
      )
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
function CliqDownMessage(csmc::CliqStateMachineContainer, status = DOWNSOLVED)

  #JT TODO maybe use Tx buffer
  newDwnMsgs = LikelihoodMessage(;
    sender = (; id = csmc.cliq.id.value, step = csmc._csm_iter),
    status = status,
  )

  # create all messages from subfg
  for mk in getCliqFrontalVarIds(csmc.cliq)
    v = getVariable(csmc.cliqSubFg, mk)
    if isInitialized(v, csmc.solveKey)
      newDwnMsgs.belief[mk] = TreeBelief(v, csmc.solveKey)
    end
  end

  logCSM(csmc, "cliq $(csmc.cliq.id), CliqDownMessage, allkeys=$(keys(newDwnMsgs.belief))")

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
  if length(getParent(csmc.tree, csmc.cliq)) != 0 &&
     getCliqueStatus(csmc.cliq) != MARGINALIZED
    logCSM(
      csmc,
      "CSM-4a doCliqDownSolve_StateMachine -- dwnmsgs=$(collect(keys(dwnmsgs.belief)))",
    )
    # maybe cycle through separators (or better yet, just use values directly -- see next line)
    msgfcts = addMsgFactors!(csmc.cliqSubFg, dwnmsgs, DownwardPass)
    # force separator variables in cliqSubFg to adopt down message values
    updateSubFgFromDownMsgs!(csmc.cliqSubFg, dwnmsgs, getCliqSeparatorVarIds(csmc.cliq))

    if dwnmsgs.status in [DOWNSOLVED, MARGINALIZED]
      logCSM(csmc, "CSM-4a doCliqDownSolve_StateMachine")
      return solveDown_StateMachine
    elseif dwnmsgs.status == INITIALIZED || dwnmsgs.status == NO_INIT
      return tryDownInit_StateMachine
    else
      logCSM(
        csmc,
        "CSM-4a Unknown state";
        status = dwnmsgs.status,
        loglevel = Logging.Error,
        c = csmc.cliqId,
      )
      error("CSM-4a waitForDown State Error: Unknown/unimplemented transision.")
    end
  else
    # Special root case or MARGINALIZED
    #TODO improve
    solveStatus = getCliqueStatus(csmc.cliq)
    logCSM(csmc, "CSM-4a root case or MARGINALIZED"; status = solveStatus, c = csmc.cliqId)
    if solveStatus in [INITIALIZED, NO_INIT, UPSOLVED, UPRECYCLED, MARGINALIZED]
      solveStatus == MARGINALIZED && setCliqueDrawColor!(csmc.cliq, "blue")
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
  dwnkeys_ =
    lsf(csmc.cliqSubFg; tags = [:__DOWNWARD_COMMON__;]) .|> x -> ls(csmc.cliqSubFg, x)[1]
  initorder = getCliqInitVarOrderDown(csmc.cliqSubFg, csmc.cliq, dwnkeys_)
  # initorder = getCliqVarInitOrderUp(csmc.tree, csmc.cliq)

  someInit = cycleInitByVarOrder!(
    csmc.cliqSubFg,
    initorder;
    solveKey = csmc.solveKey,
    logger = csmc.logger,
  )
  # is clique fully upsolved or only partially?
  # print out the partial init status of all vars in clique
  printCliqInitPartialInfo(csmc.cliqSubFg, csmc.cliq, csmc.solveKey, csmc.logger)
  logCSM(csmc, "CSM-4b tryInitCliq_StateMachine -- someInit=$someInit, varorder=$initorder")

  msgfcts = deleteMsgFactors!(csmc.cliqSubFg, [:__DOWNWARD_COMMON__;]) # msgfcts # TODO, use tags=[:__LIKELIHOODMESSAGE__], see #760
  logCSM(
    csmc,
    "CSM-4b tryDownInit_StateMachine - removing factors, length=$(length(msgfcts))",
  )

  solveStatus = someInit ? INITIALIZED : NO_INIT
  if someInit
    setCliqueDrawColor!(csmc.cliq, "seagreen")
  else
    setCliqueDrawColor!(csmc.cliq, "khaki")
  end
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

  # get down msg from Rx buffer (saved in take!)
  dwnmsgs = getMessageBuffer(csmc.cliq).downRx
  # #XXX
  # # maybe cycle through separators (or better yet, just use values directly -- see next line)
  # msgfcts = addMsgFactors!(csmc.cliqSubFg, dwnmsgs, DownwardPass)  
  # force separator variables in cliqSubFg to adopt down message values
  # updateSubFgFromDownMsgs!(csmc.cliqSubFg, dwnmsgs, getCliqSeparatorVarIds(csmc.cliq))

  opts = getSolverParams(csmc.dfg)
  #XXX test with and without
  # add required all frontal connected factors
  if !opts.useMsgLikelihoods
    newvars, newfcts = addDownVariableFactors!(
      csmc.dfg,
      csmc.cliqSubFg,
      csmc.cliq,
      csmc.logger;
      solvable = 1,
    )
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

  logCSM(csmc, "CSM-4c $(csmc.cliq.id): clique down solve completed")

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
    logCSM(
      csmc,
      "CSM-4d Empty message on clique $(csmc.cliqId) frontals";
      loglevel = Logging.Info,
    )
  end

  logCSM(
    csmc,
    "CSM-4d msg to send down on $(keys(beliefMsg.belief))";
    beliefMsg = beliefMsg,
  )
  # pass through the frontal variables that were sent from above
  downmsg = getMessageBuffer(csmc.cliq).downRx
  svars = getCliqSeparatorVarIds(csmc.cliq)
  if !isnothing(downmsg)
    pass_through_separators = intersect(svars, keys(downmsg.belief))
    for si in pass_through_separators
      beliefMsg.belief[si] = downmsg.belief[si]
      logCSM(csmc, "CSM-4d adding parent message"; sym = si, msg = downmsg.belief[si])
    end
  end

  # Store the down message for debugging, will be stored even if no children present
  getMessageBuffer(csmc.cliq).downTx = beliefMsg

  #TODO maybe send a specific message to only the child that needs it
  @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
    logCSM(csmc, "CSM-4d $(csmc.cliq.id): put! on edge $(e)")
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
  isParametricSolve = csmc.algorithm == :parametric

  # set PPE and solved for all frontals
  if !isParametricSolve
    for sym in getCliqFrontalVarIds(csmc.cliq)
      # set PPE in cliqSubFg
      setVariablePosteriorEstimates!(csmc.cliqSubFg, sym)
      # set solved flag
      vari = getVariable(csmc.cliqSubFg, sym, csmc.solveKey)
      setSolvedCount!(vari, getSolvedCount(vari, csmc.solveKey) + 1, csmc.solveKey)
    end
  end

  # transfer results to main factor graph
  frsyms = getCliqFrontalVarIds(csmc.cliq)
  logCSM(
    csmc,
    "CSM-5 finishingCliq -- transferUpdateSubGraph! with solveKey=$(csmc.solveKey) on $frsyms",
  )
  transferUpdateSubGraph!(
    csmc.dfg,
    csmc.cliqSubFg,
    frsyms,
    csmc.logger;
    solveKey = csmc.solveKey,
    updatePPE = !isParametricSolve,
  )

  #solve finished change color
  setCliqueDrawColor!(csmc.cliq, "lightblue")

  logCSM(
    csmc,
    "CSM-5 Clique $(csmc.cliq.id) finished, solveKey=$(csmc.solveKey)";
    loglevel = Logging.Debug,
  )
  return IncrementalInference.exitStateMachine
end

#
