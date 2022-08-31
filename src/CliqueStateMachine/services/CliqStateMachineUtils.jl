
## ===================================================================================================================
##  Consolidate by first changing to common functions of csmc
## ===================================================================================================================
"""
    $SIGNATURES

Calculate the full upward Chapman-Kolmogorov transit integral solution approximation (i.e. upsolve).

Notes
- State machine function nr. 8g
- Assumes LIKELIHOODMESSAGE factors are in csmc.cliqSubFg but does not remove them.
- TODO: Make multi-core

DevNotes
- NEEDS DFG v0.8.1, see IIF #760
- temperory consolidation function
"""
function __doCliqUpSolveInitialized!(csmc::CliqStateMachineContainer)

  # check if all cliq vars have been initialized so that full inference can occur on clique
  status = getCliqueStatus(csmc.cliq)
  infocsm(csmc, "8g, doCliqUpSolveInitialized_StateMachine -- clique status = $(status)")
  logCSM(csmc, "8g, doCliqUpSolveInitialized_StateMachine -- clique status = $(status)")

  setCliqueDrawColor!(csmc.cliq, "red")

  opt = getSolverParams(csmc.cliqSubFg)
  # get Dict{Symbol, TreeBelief} of all updated variables in csmc.cliqSubFg
  retdict = approxCliqMarginalUp!(csmc; iters = opt.gibbsIters, logger = csmc.logger)
  # retdict = approxCliqMarginalUp!(csmc, LikelihoodMessage[]; iters=4, logger=csmc.logger)
  logCSM(csmc, "aproxCliqMarginalUp!"; retdict = retdict)

  updateFGBT!(
    csmc.cliqSubFg,
    csmc.cliq,
    retdict;
    dbg = getSolverParams(csmc.cliqSubFg).dbg,
    logger = csmc.logger,
  ) # urt

  # set clique color accordingly, using local memory
  # setCliqueDrawColor!(csmc.cliq, isCliqFullDim(csmc.cliqSubFg, csmc.cliq) ? "pink" : "tomato1")

  # notify of results (part of #459 consolidation effort)
  getCliqueData(csmc.cliq).upsolved = true

  return nothing
end

## ===================================================================================================================
##  CSM logging functions
## ===================================================================================================================
# using Serialization

"""
    $SIGNATURES

Internal helper function to save a dfg object to LogPath during clique state machine operations.

Notes
- will only save dfg object if `opts.dbg=true`

Related

saveDFG, loadDFG!, loadDFG
"""
function _dbgCSMSaveSubFG(csmc::CliqStateMachineContainer, filename::String)
  opt = getSolverParams(csmc.cliqSubFg)

  if opt.dbg
    folder::String = joinpath(opt.logpath, "logs", "cliq$(getId(csmc.cliq))")
    if !ispath(folder)
      mkpath(folder)
    end
    # NOTE there was a bug using saveDFG, so used serialize, left for future use  
    # serialize(joinpath(folder, filename), csmc.cliqSubFg)
    DFG.saveDFG(csmc.cliqSubFg, joinpath(folder, filename))
    drawGraph(csmc.cliqSubFg; show = false, filepath = joinpath(folder, "$(filename).pdf"))
  end

  return opt.dbg
end

"""
    $SIGNATURES

Specialized info logger print function to show clique state machine information
in a standardized form.
"""
function infocsm(csmc::CliqStateMachineContainer, str::A) where {A <: AbstractString}
  tm = string(Dates.now())
  tmt = split(tm, 'T')[end]

  lbl = getLabel(csmc.cliq)
  lbl1 = split(lbl, ',')[1]
  cliqst = getCliqueStatus(csmc.cliq)

  with_logger(csmc.logger) do
    @info "$tmt | $(getId(csmc.cliq))---$lbl1 @ $(cliqst) | " * str
  end
  flush(csmc.logger.stream)
  return nothing
end

"""
    $SIGNATURES
Helper function to log a message at a specific level to a clique identified by `csm_i` where i = cliq.id
Notes:
- Related to infocsm.
- Different approach to logging that uses the build in logging functionality to provide more flexibility.
- Can be used with LoggingExtras.jl 
"""
function logCSM(
  csmc,
  msg::String;
  loglevel::Logging.LogLevel = Logging.Debug,
  maxlog = nothing,
  kwargs...,
)
  csmc.enableLogging ? nothing : (return nothing)
  #Debug = -1000
  #Info = 0
  #Warn = 1000
  #Error = 2000
  @logmsg(
    loglevel,
    msg,
    _module = begin
      bt = backtrace()
      funcsym = (:logCSM, Symbol("logCSM##kw")) #always use the calling function of logCSM
      frame, caller = Base.firstcaller(bt, funcsym)
      # TODO: Is it reasonable to attribute callers without linfo to Core?
      caller.linfo isa Core.MethodInstance ? caller.linfo.def.module : Core
    end,
    _file = String(caller.file),
    _line = caller.line,
    _id = (frame, funcsym),
    # caller=caller, 
    # st4 = stacktrace()[4],
    _group = Symbol("csm_$(csmc.cliq.id)"),
    maxlog = maxlog,
    kwargs...
  )

  return nothing
end

## ===================================================================================================================
## CSM Error functions
## ===================================================================================================================

function putErrorDown(csmc::CliqStateMachineContainer)
  setCliqueDrawColor!(csmc.cliq, "red")
  @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
    logCSM(csmc, "CSM clique $(csmc.cliq.id): propagate down error on edge $(e)")
    @async putBeliefMessageDown!(csmc.tree, e, LikelihoodMessage(; status = ERROR_STATUS))
  end
  logCSM(
    csmc,
    "CSM clique $(csmc.cliq.id): Exit with error state";
    loglevel = Logging.Error,
  )
  return nothing
end

function putErrorUp(csmc::CliqStateMachineContainer)
  setCliqueDrawColor!(csmc.cliq, "red")
  for e in getEdgesParent(csmc.tree, csmc.cliq)
    logCSM(csmc, "CSM clique, $(csmc.cliq.id): propagate up error on edge $(e)")
    putBeliefMessageUp!(csmc.tree, e, LikelihoodMessage(; status = ERROR_STATUS))
  end
  return nothing
end

## ===================================================================================================================
## CSM Monitor functions
## ===================================================================================================================

"""
    $SIGNATURES
Monitor CSM tasks for failures and propagate error to the other CSMs to cleanly exit. 
"""
function monitorCSMs(tree, alltasks; forceIntExc::Bool = false)
  task = @async begin
    while true
      all(istaskdone.(alltasks)) && (@info "monitorCSMs: all tasks done"; break)
      for (i, t) in enumerate(alltasks)
        if istaskfailed(t)
          if forceIntExc
            @error "Task $i failed, sending InterruptExceptions to all running CSM tasks"
            throwIntExcToAllTasks(alltasks)
            @debug "done with throwIntExcToAllTasks"
          else
            @error "Task $i failed, sending error to all cliques"
            bruteForcePushErrorCSM(tree)
            # for tree.messageChannels
            @info "All cliques should have exited"
          end
        end
      end
      sleep(1)
    end
  end
  return task
end

function throwIntExcToAllTasks(alltasks)
  for (i, t) in enumerate(alltasks)
    if !istaskdone(alltasks[i])
      @debug "Sending InterruptExceptions to CSM task $i"
      schedule(alltasks[i], InterruptException(); error = true)
      @debug "InterruptExceptions CSM task $i"
    end
  end
  return nothing
end

function bruteForcePushErrorCSM(tree::AbstractBayesTree)
  errMsg = LikelihoodMessage(; status = ERROR_STATUS)
  for (i, ch) in getMessageChannels(tree)
    if isready(ch.upMsg)
      take!(ch.upMsg)
    else
      @debug("Up edge $i", ch.upMsg)
      @async put!(ch.upMsg, errMsg)
    end
    if isready(ch.downMsg)
      take!(ch.downMsg)
    else
      @debug("Down edge $i", ch.downMsg)
      @async put!(ch.downMsg, errMsg)
    end
  end

  for (i, ch) in getMessageChannels(tree)
    while isready(ch.upMsg)
      @debug "cleanup take on $i up"
      take!(ch.upMsg)
    end
    while isready(ch.downMsg)
      @debug "cleanup take on $i down"
      take!(ch.downMsg)
    end
  end
end

## ===================================================================================================================
## CSM Clique Functions
## ===================================================================================================================

"""
    $SIGNATURES

Set all up `upsolved` and `downsolved` cliq data flags `to::Bool=false`.
"""
function setAllSolveFlags!(treel::AbstractBayesTree, to::Bool = false)::Nothing
  for (id, cliq) in getCliques(treel)
    cliqdata = getCliqueData(cliq)
    setCliqueStatus!(cliqdata, NULL)
    cliqdata.upsolved = to
    cliqdata.downsolved = to
  end
  return nothing
end

"""
    $SIGNATURES

Return true or false depending on whether the tree has been fully initialized/solved/marginalized.
"""
function isTreeSolved(treel::AbstractBayesTree; skipinitialized::Bool = false)
  acclist = CliqStatus[UPSOLVED, DOWNSOLVED, MARGINALIZED]
  skipinitialized ? nothing : push!(acclist, INITIALIZED)
  for (clid, cliq) in getCliques(treel)
    if !(getCliqueStatus(cliq) in acclist)
      return false
    end
  end
  return true
end

function isTreeSolvedUp(treel::AbstractBayesTree)
  for (clid, cliq) in getCliques(treel)
    if getCliqueStatus(cliq) != UPSOLVED
      return false
    end
  end
  return true
end

"""
    $SIGNATURES

Reset the Bayes (Junction) tree so that a new upsolve can be performed.

Notes
- Will change previous clique status from `DOWNSOLVED` to `INITIALIZED` only.
- Sets the color of tree clique to `lightgreen`.
"""
function resetTreeCliquesForUpSolve!(treel::AbstractBayesTree)::Nothing
  acclist = CliqStatus[DOWNSOLVED]
  for (clid, cliq) in getCliques(treel)
    if getCliqueStatus(cliq) in acclist
      setCliqueStatus!(cliq, INITIALIZED)
      setCliqueDrawColor!(cliq, "sienna")
    end
  end
  return nothing
end

"""
    $SIGNATURES

Return true there is no other sibling that will make progress.

Notes
- Relies on sibling priority order with only one "currently best" option that will force progress in global upward inference.
- Return false if one of the siblings is still busy
"""
function areSiblingsRemaingNeedDownOnly(tree::AbstractBayesTree, cliq::TreeClique)::Bool
  #
  stillbusylist = [NULL, INITIALIZED]
  prnt = getParent(tree, cliq)
  if length(prnt) > 0
    for si in getChildren(tree, prnt[1])
      # are any of the other siblings still busy?
      if si.id != cliq.id && getCliqueStatus(si) in stillbusylist
        return false
      end
    end
  end

  # nope, everybody is waiting for something to change -- proceed with forcing a cliq solve
  return true
end

"""
    $SIGNATURES

Approximate Chapman-Kolmogorov transit integral and return separator marginals as messages to pass up the Bayes (Junction) tree, along with additional clique operation values for debugging.

Notes
- `onduplicate=true` by default internally uses deepcopy of factor graph and Bayes tree, and does **not** update the given objects.  Set false to update `fgl` and `treel` during compute.

Future
- TODO: internal function chain is too long and needs to be refactored for maintainability.
"""
function approxCliqMarginalUp!(
  csmc::CliqStateMachineContainer,
  childmsgs = LikelihoodMessage[];#fetchMsgsUpChildren(csmc, TreeBelief);
  N::Int = getSolverParams(csmc.cliqSubFg).N,
  dbg::Bool = getSolverParams(csmc.cliqSubFg).dbg,
  multiproc::Bool = getSolverParams(csmc.cliqSubFg).multiproc,
  logger = ConsoleLogger(),
  iters::Int = 3,
  drawpdf::Bool = false,
)
  #
  # use subgraph copy of factor graph for operations and transfer variables results later only
  fg_ = csmc.cliqSubFg
  tree_ = csmc.tree
  cliq = csmc.cliq

  with_logger(logger) do
    @info "======== Clique $(getLabel(cliq)) ========"
  end

  if multiproc
    cliqc = deepcopy(cliq)
    btnd = getCliqueData(cliqc)
    # ett.cliq = cliqc
    # TODO create new dedicated file for separate process to log with
    try
      retdict = remotecall_fetch(
        upGibbsCliqueDensity,
        getWorkerPool(),
        fg_,
        cliqc,
        csmc.solveKey,
        childmsgs,
        N,
        dbg,
        iters,
      )
    catch ex
      with_logger(logger) do
        @info ex
        @error ex
        flush(logger.stream)
        msg = sprint(showerror, ex)
        @error msg
      end
      flush(logger.stream)
      error(ex)
    end
  else
    with_logger(logger) do
      @info "Single process upsolve clique=$(cliq.id)"
    end
    retdict =
      upGibbsCliqueDensity(fg_, cliq, csmc.solveKey, childmsgs, N, dbg, iters, logger)
  end

  with_logger(logger) do
    @info "=== end Clique $(getLabel(cliq)) ========================"
  end

  return retdict
end

"""
    $SIGNATURES

Determine which variables to iterate or compute directly for downward tree pass of inference.

DevNotes
- # TODO see #925

Related

directPriorMsgIDs, directFrtlMsgIDs, directAssignmentIDs, mcmcIterationIDs
"""
function determineCliqVariableDownSequence(
  subfg::AbstractDFG,
  cliq::TreeClique;
  solvable::Int = 1,
  logger = ConsoleLogger(),
)
  #
  frtl = getCliqFrontalVarIds(cliq)

  adj, varLabels, FactorLabels = DFG.getBiadjacencyMatrix(subfg; solvable = solvable)
  mask = (x -> x in frtl).(varLabels)
  newFrtlOrder = varLabels[mask]
  subAdj = adj[:, mask]
  #TODO don't use this getAdjacencyMatrixSymbols, #604
  # adj = DFG.getAdjacencyMatrixSymbols(subfg, solvable=solvable)
  # mask = map(x->(x in frtl), adj[1,:])
  # subAdj = adj[2:end,mask] .!= nothing
  # newFrtlOrder = Symbol.(adj[1,mask])

  crossCheck = 1 .< sum(Int.(subAdj); dims = 2)
  iterVars = Symbol[]
  for i = 1:length(crossCheck)
    # must add associated variables to iterVars
    if crossCheck[i]
      # # DEBUG loggin
      # with_logger(logger) do
      #   @info "newFrtlOrder=$newFrtlOrder"
      #   @info "(subAdj[i,:]).nzind=$((subAdj[i,:]).nzind)"
      # end
      # flush(logger.stream)
      # find which variables are associated
      varSym = newFrtlOrder[(subAdj[i, :]).nzind]
      union!(iterVars, varSym)
    end
  end

  # return iteration list ordered by frtl
  return intersect(frtl, iterVars)
end

"""
    $SIGNATURES

Perform downward direction solves on a sub graph fragment.
Calculates belief on each of the frontal variables and iterate if required.

Notes
- uses all factors connected to the frontal variables.
- assumes `subfg` was properly prepared before calling.
- has multi-process option.

Dev Notes
- TODO incorporate variation possible due to cross frontal factors.
- cleanup and updates required, and @spawn jl 1.3
"""
function solveCliqDownFrontalProducts!(
  subfg::AbstractDFG,
  cliq::TreeClique,
  opts::SolverParams,
  logger = ConsoleLogger();
  solveKey::Symbol = :default,
  MCIters::Int = 3,
)
  #
  # get frontal variables for this clique
  frsyms = getCliqFrontalVarIds(cliq)

  # determine if cliq has cross frontal factors
  # iterdwn, directdwns, passmsgs?
  iterFrtls = determineCliqVariableDownSequence(subfg, cliq; logger = logger)

  # direct frontals
  directs = setdiff(frsyms, iterFrtls)

  # ignore limited fixed lag variables
  fixd = map(x -> opts.limitfixeddown && isMarginalized(subfg, x), frsyms)
  skip = frsyms[fixd]
  iterFrtls = setdiff(iterFrtls, skip)
  directs = setdiff(directs, skip)
  with_logger(logger) do
    @info "cliq $(cliq.id), solveCliqDownFrontalProducts!, skipping marginalized keys=$(skip)"
  end

  # use new localproduct approach
  if opts.multiproc
    downresult =
      Dict{Symbol, Tuple{ManifoldKernelDensity, Vector{Float64}, Vector{Symbol}}}()
    @sync for i = 1:length(directs)
      @async begin
        downresult[directs[i]] = remotecall_fetch(
          localProductAndUpdate!,
          getWorkerPool(),
          subfg,
          directs[i],
          false;
          solveKey = solveKey,
        )
        # downresult[directs[i]] = remotecall_fetch(localProductAndUpdate!, upp2(), subfg, directs[i], false)
      end
    end
    with_logger(logger) do
      @info "cliq $(cliq.id), solveCliqDownFrontalProducts!, multiproc keys=$(keys(downresult))"
    end
    for fr in directs
      with_logger(logger) do
        @info "cliq $(cliq.id), solveCliqDownFrontalProducts!, key=$(fr), infdim=$(downresult[fr][2]), lbls=$(downresult[fr][3])"
      end
      setValKDE!(subfg, fr, downresult[fr][1], false, downresult[fr][2])
    end
    for mc = 1:MCIters, fr in iterFrtls
      try
        result = remotecall_fetch(
          localProductAndUpdate!,
          getWorkerPool(),
          subfg,
          fr,
          false;
          solveKey = solveKey,
        )
        # result = remotecall_fetch(localProductAndUpdate!, upp2(), subfg, fr, false)
        setValKDE!(subfg, fr, result[1], false, result[2])
        with_logger(logger) do
          @info "cliq $(cliq.id), solveCliqDownFrontalProducts!, iter key=$(fr), infdim=$(result[2]), lbls=$(result[3])"
        end
      catch ex
        # what if results contains an error?
        with_logger(logger) do
          @error ex
          flush(logger.stream)
          msg = sprint(showerror, ex)
          @error msg
        end
        error(ex)
      end
    end
  else
    # do directs first
    for fr in directs
      localProductAndUpdate!(subfg, fr, true, logger; solveKey = solveKey)
    end
    #do iters next
    for mc = 1:MCIters, fr in iterFrtls
      localProductAndUpdate!(subfg, fr, true, logger; solveKey = solveKey)
    end
  end

  return nothing
end

#
