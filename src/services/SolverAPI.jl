## Various solver API's used in the past.  These functions are due to be standardized, and obsolete code / functions removed.

export solveTree!, solveGraph!
export fetchCliqHistoryAll!

## ==============================================================================================
## Launch the tasks/treads for cliques
## ==============================================================================================
"""
    $SIGNATURES

Start tasks (@async or Threads.@spawn threads if multithread=true) to solve the factor graph on the tree.
"""
function taskSolveTree!(
  dfg::AbstractDFG,
  treel::AbstractBayesTree,
  timeout::Union{Nothing, <:Real} = nothing;
  oldtree::AbstractBayesTree = BayesTree(),
  drawtree::Bool = false,
  verbose::Bool = false,
  verbosefid = stdout,
  limititers::Int = -1,
  limititercliqs::Vector{Pair{Symbol, Int}} = Pair{Symbol, Int}[],
  downsolve::Bool = false,
  incremental::Bool = false,
  multithread::Bool = false,
  skipcliqids::Vector{Symbol} = Symbol[],
  recordcliqs::Vector{Symbol} = Symbol[],
  delaycliqs::Vector{Symbol} = Symbol[],
  smtasks = Task[],
  algorithm::Symbol = :default,
  solveKey::Symbol = algorithm,
)
  #
  # revert DOWNSOLVED status to INITIALIZED in preparation for new upsolve
  resetTreeCliquesForUpSolve!(treel)

  drawtree ? drawTree(treel; show = false, filepath = joinLogPath(dfg, "bt.dot")) : nothing

  cliqHistories = Dict{Int, Vector{CSMHistoryTuple}}()

  resize!(smtasks, getNumCliqs(treel))

  approx_iters = getNumCliqs(treel) * 24
  solve_progressbar =
    verbose ? nothing : ProgressUnknown("Solve Progress: approx max $approx_iters, at iter")

  # queue all the tasks/threads
  if !isTreeSolved(treel; skipinitialized = true)
    @sync begin
      monitortask = monitorCSMs(treel, smtasks)
      # duplicate int i into async (important for concurrency)
      for i = 1:getNumCliqs(treel) # TODO, this might not always work?
        scsym = getCliqFrontalVarIds(getClique(treel, i))
        if length(intersect(scsym, skipcliqids)) == 0
          limthiscsm = filter(x -> (x[1] in scsym), limititercliqs)
          limiter = 0 < length(limthiscsm) ? limthiscsm[1][2] : limititers

          if multithread
            smtasks[i] = Threads.@spawn tryCliqStateMachineSolve!(
              dfg,
              treel,
              i,
              timeout;
              solveKey = solveKey,
              algorithm = algorithm,
              oldtree = oldtree,
              verbose = verbose,
              verbosefid = verbosefid,
              drawtree = drawtree,
              limititers = limititers,
              downsolve = downsolve,
              incremental = incremental,
              delaycliqs = delaycliqs,
              recordcliqs = recordcliqs,
              solve_progressbar = solve_progressbar,
            )
          else
            smtasks[i] = @async tryCliqStateMachineSolve!(
              dfg,
              treel,
              i,
              timeout;
              solveKey = solveKey,
              algorithm = algorithm,
              oldtree = oldtree,
              verbose = verbose,
              verbosefid = verbosefid,
              drawtree = drawtree,
              limititers = limiter,
              downsolve = downsolve,
              incremental = incremental,
              delaycliqs = delaycliqs,
              recordcliqs = recordcliqs,
              solve_progressbar = solve_progressbar,
            )
          end
        end # if
      end # for
    end # sync
  end # if

  # if record cliques is in use, else skip computational delay
  0 == length(recordcliqs) ? nothing : fetchCliqHistoryAll!(smtasks, cliqHistories)

  !isnothing(solve_progressbar) && finish!(solve_progressbar)

  return smtasks, cliqHistories
end

function tryCliqStateMachineSolve!(
  dfg::G,
  treel::AbstractBayesTree,
  cliqKey::Union{Int, CliqueId},
  timeout::Union{Nothing, <:Real} = nothing;
  oldtree::AbstractBayesTree = BayesTree(),
  verbose::Bool = false,
  verbosefid = stdout,
  drawtree::Bool = false,
  limititers::Int = -1,
  downsolve::Bool = false,
  incremental::Bool = false,
  delaycliqs::Vector{Symbol} = Symbol[],
  recordcliqs::Vector{Symbol} = Symbol[],
  solve_progressbar = nothing,
  algorithm::Symbol = :default,
  solveKey::Symbol = algorithm,
) where {G <: AbstractDFG}
  #
  clst = :na
  cliq = getClique(treel, cliqKey)
  syms = getCliqFrontalVarIds(cliq)

  oldcliq = attemptTreeSimilarClique(oldtree, getCliqueData(cliq))
  oldcliqdata = getCliqueData(oldcliq)

  opts = getSolverParams(dfg)
  # Base.rm(joinpath(opts.logpath,"logs/cliq$i"), recursive=true, force=true)
  mkpath(joinpath(opts.logpath, "logs/cliq$(cliq.id)/"))
  logger = SimpleLogger(open(joinpath(opts.logpath, "logs/cliq$(cliq.id)/log.txt"), "w+")) # NullLogger()
  # global_logger(logger)
  history = Vector{CSMHistoryTuple}()
  recordthiscliq = length(intersect(recordcliqs, syms)) > 0
  delaythiscliq = length(intersect(delaycliqs, syms)) > 0
  try
    history = initStartCliqStateMachine!(
      dfg,
      treel,
      cliq,
      timeout;
      oldcliqdata = oldcliqdata,
      drawtree = drawtree,
      verbose = verbose,
      verbosefid = verbosefid,
      limititers = limititers,
      downsolve = downsolve,
      recordhistory = recordthiscliq,
      incremental = incremental,
      delay = delaythiscliq,
      logger = logger,
      solve_progressbar = solve_progressbar,
      algorithm = algorithm,
      solveKey = solveKey,
    )
    #
    # cliqHistories[cliqKey] = history
    if length(history) >= limititers && limititers != -1
      # @warn "writing logs/cliq$(cliq.id)/csm.txt"
      # @save "/tmp/cliqHistories/cliq$(cliq.id).jld2" history
      fid = open(joinpath(opts.logpath, "logs/cliq$(cliq.id)/csm.txt"), "w")
      printCliqHistorySummary(fid, history)
      close(fid)
    end
    flush(logger.stream)
    close(logger.stream)
    # clst = getCliqueStatus(cliq)
    # clst = cliqInitSolveUp!(dfg, treel, cliq, drawtree=drawtree, limititers=limititers )
  catch err
    bt = catch_backtrace()
    println()
    showerror(stderr, err, bt)
    # @warn "writing /tmp/caesar/logs/cliq$(cliq.id)/*.txt"
    fid = open(joinpath(opts.logpath, "logs/cliq$(cliq.id)/stacktrace.txt"), "w")
    showerror(fid, err, bt)
    close(fid)
    fid = open(joinpath(opts.logpath, "logs/cliq$(cliq.id)_stacktrace.txt"), "w")
    showerror(fid, err, bt)
    close(fid)
    # @save "/tmp/cliqHistories/$(cliq.label).jld2" history
    fid = open(joinpath(opts.logpath, "logs/cliq$(cliq.id)/csm.txt"), "w")
    printCliqHistorySummary(fid, history)
    close(fid)
    fid = open(joinpath(opts.logpath, "logs/cliq$(cliq.id)_csm.txt"), "w")
    printCliqHistorySummary(fid, history)
    close(fid)
    flush(logger.stream)
    close(logger.stream)
    rethrow()
  end
  # if !(clst in [UPSOLVED; DOWNSOLVED; MARGINALIZED])
  #   error("Clique $(cliq.id), initInferTreeUp! -- cliqInitSolveUp! did not arrive at the desired solution statu: $clst")
  # end
  return history
end

"""
    $SIGNATURES

Standalone state machine solution for a single clique.

Related:

initInferTreeUp!
"""
function solveCliqWithStateMachine!(
  dfg::G,
  tree::AbstractBayesTree,
  frontal::Symbol;
  iters::Int = 200,
  downsolve::Bool = true,
  recordhistory::Bool = false,
  verbose::Bool = false,
  nextfnc::Function = canCliqMargRecycle_StateMachine,
  prevcsmc::Union{Nothing, CliqStateMachineContainer} = nothing,
) where {G <: AbstractDFG}
  #
  cliq = getClique(tree, frontal)

  children = getChildren(tree, cliq)#Graphs.out_neighbors(cliq, tree.bt)

  prnt = getParent(tree, cliq)

  destType = (G <: InMemoryDFGTypes) ? G : LocalDFG

  csmc = if isa(prevcsmc, Nothing)
    CliqStateMachineContainer(
    dfg,
    initfg(destType; solverParams = getSolverParams(dfg)),
    tree,
    cliq,
    prnt,
    children,
    false,
    true,
    true,
    downsolve,
    false,
    getSolverParams(dfg),
  )
  else
    prevcsmc
  end
  statemachine =
    StateMachine{CliqStateMachineContainer}(; next = nextfnc, name = "cliq$(cliq.id)")
  while statemachine(
    csmc;
    verbose = verbose,
    iterlimit = iters,
    recordhistory = recordhistory,
  )
  end
  return statemachine, csmc
end

## ==============================================================================================
# Prepare CSM (based on FSM) entry points
## ==============================================================================================

"""
    $SIGNATURES

Fetch solver history from clique state machines that have completed their async Tasks and store in the `hist::Dict{Int,Tuple}` dictionary.
"""
function fetchCliqHistoryAll!(
  smt::Vector{Task},
  hist::Dict{Int, Vector{CSMHistoryTuple}} = Dict{Int, Vector{CSMHistoryTuple}}(),
)
  #
  for i = 1:length(smt)
    sm = smt[i]
    # only fetch states that have completed processing
    if sm.state == :done
      haskey(hist, i) ? @warn("overwriting existing history key $i") : nothing
      hist[i] = fetch(sm)
    elseif !isnothing(sm.storage) && haskey(sm.storage, :statemachine)
      hist[i] = CSMHistoryTuple.(sm.storage[:statemachine].history)
    end
  end
  return hist
end

## ==============================================================================================
# Nominal user interface to the solver
## ==============================================================================================

"""
    $SIGNATURES

Perform inference over the Bayes tree according to `opt::SolverParams` and keyword arguments.

Notes
- Aliased with `solveGraph!`
- Variety of options, including fixed-lag solving -- see `getSolverParams(fg)` for details.
  - See online Documentation for more details: https://juliarobotics.org/Caesar.jl/latest/
- Latest result always stored in `solvekey=:default`.
- Experimental `storeOld::Bool=true` will duplicate the current result as supersolve `:default_k`.
  - Based on `solvable==1` assumption.
- `limititercliqs` allows user to limit the number of iterations a specific CSM does.
- keywords `verbose` and `verbosefid::IOStream` can be used together to to send output to file or default `stdout`.
- keyword `recordcliqs=[:x0; :x7...]` identifies by frontals which cliques to record CSM steps.
  - See [`repeatCSMStep!`](@ref), [`printCSMHistoryLogical`](@ref), [`printCSMHistorySequential`](@ref)

DevNotes
- TODO Change keyword arguments to new @parameter `SolverOptions` type.

Example
```julia
# pass in old `tree` to enable compute recycling -- see online Documentation for more details
tree = solveTree!(fg [,tree])
```

Related

`solveGraph!`, [`solveCliqUp!`](@ref), [`solveCliqDown!`](@ref), [`buildTreeReset!`](@ref), [`repeatCSMStep`](@ref), [`printCSMHistoryLogical`](@ref)
"""
function solveTree!(
  dfgl::AbstractDFG,
  oldtree::AbstractBayesTree = BayesTree();
  timeout::Union{Nothing, <:Real} = nothing,
  storeOld::Bool = false,
  verbose::Bool = false,
  verbosefid = stdout,
  delaycliqs::Vector{Symbol} = Symbol[],
  recordcliqs::Vector{Symbol} = Symbol[],
  limititercliqs::Vector{Pair{Symbol, Int}} = Pair{Symbol, Int}[],
  injectDelayBefore::Union{Nothing, Vector{<:Pair{Int, <:Pair{<:Function, <:Real}}}} = nothing,
  skipcliqids::Vector{Symbol} = Symbol[],
  eliminationOrder::Union{Nothing, Vector{Symbol}} = nothing,
  eliminationConstraints::Vector{Symbol} = Symbol[],
  smtasks::Vector{Task} = Task[],
  dotreedraw = Int[1;],
  runtaskmonitor::Bool = true,
  algorithm::Symbol = :default,
  solveKey::Symbol = algorithm,
  multithread::Bool = true,
)
  #
  # workaround in case isolated variables occur
  ensureSolvable!(dfgl)
  opt = getSolverParams(dfgl)

  # showtree should force drawtree
  if opt.showtree && !opt.drawtree
    @info("Since .showtree=true, also bumping .drawtree=true")
  else
    nothing
  end
  opt.drawtree |= opt.showtree

  # depcrecation
  # update worker pool incase there are more or less
  setWorkerPool!()
  if opt.multiproc && nprocs() == 1
    @info "Setting `.multiproc=false` since `Distributed.nprocs() == 1`"
    opt.multiproc = false
  end
  
  if opt.graphinit
    @info "Ensure variables are all initialized (graphinit)"
    if algorithm == :parametric
      @warn "Parametric is using default graphinit (and ignoring solveKey)"
      initAll!(dfgl)
      initParametricFrom!(dfgl)
    else
      initAll!(dfgl, solveKey)
    end
  end
  # construct tree
  @info "Solving over the Bayes (Junction) tree."

  hist = Dict{Int, Vector{CSMHistoryTuple}}()

  if opt.isfixedlag
    @info "Quasi fixed-lag is enabled (a feature currently in testing, and ignoring solveKey)!"
    fifoFreeze!(dfgl)
  end

  # perhaps duplicate current value
  if storeOld || opt.dbg
    ss = listSupersolves(dfgl) .|> string
    ss_ = ss[occursin.(r"default_", ss)] .|> x -> x[9:end]
    filter!(x -> occursin(r"^\d+$", x), ss_)  # ss_ = ss_[occursin.(r"^\d$",ss_)]
    allk = parse.(Int, ss_)
    nextk = length(allk) == 0 ? 0 : maximum(allk) + 1
    newKey = Symbol(:default_, nextk)
    cloneSolveKey!(dfgl, newKey, :default; solvable = 1)
    # foreach(x->updateVariableSolverData!(dfgl, x, getSolverData(getVariable(dfgl,x), :default), newKey, true, Symbol[]), ls(dfgl, solvable=1))
    @info "storeOld=true, previous :default deepcopied into $newKey for solvable==1 variables."
  end

  orderMethod = 0 < length(eliminationConstraints) ? :ccolamd : :qr

  # current incremental solver builds a new tree and matches against old tree for recycling.
  tree = buildTreeReset!(
    dfgl,
    eliminationOrder;
    drawpdf = false,
    show = opt.showtree,
    ensureSolvable = false,
    filepath = joinpath(opt.logpath, "bt.pdf"),
    eliminationConstraints = eliminationConstraints,
    ordering = orderMethod,
  )

  # setAllSolveFlags!(tree, false)

  initTreeMessageChannels!(tree)

  # if desired, drawtree in a loop
  treetask, _dotreedraw = drawTreeAsyncLoop(tree, opt; dotreedraw = dotreedraw)

  @info "Do tree based init-ference"
  algorithm != :parametric ? nothing : @error("Under development, do not use, see #539")
  !storeOld ? nothing : @error("parametric storeOld keyword not wired up yet.")

  if opt.async
    @async smtasks, hist = taskSolveTree!(
      dfgl,
      tree,
      timeout;
      solveKey = solveKey,
      algorithm = algorithm,
      multithread = multithread,
      smtasks = smtasks,
      oldtree = oldtree,
      verbose = verbose,
      verbosefid = verbosefid,
      drawtree = opt.drawtree,
      recordcliqs = recordcliqs,
      limititers = opt.limititers,
      downsolve = opt.downsolve,
      incremental = opt.incremental,
      skipcliqids = skipcliqids,
      delaycliqs = delaycliqs,
      limititercliqs = limititercliqs,
    )
  else
    smtasks, hist = taskSolveTree!(
      dfgl,
      tree,
      timeout;
      solveKey = solveKey,
      algorithm = algorithm,
      multithread = multithread,
      smtasks = smtasks,
      oldtree = oldtree,
      verbose = verbose,
      verbosefid = verbosefid,
      drawtree = opt.drawtree,
      recordcliqs = recordcliqs,
      limititers = opt.limititers,
      downsolve = opt.downsolve,
      incremental = opt.incremental,
      skipcliqids = skipcliqids,
      delaycliqs = delaycliqs,
      limititercliqs = limititercliqs,
    )
    @info "Finished tree based init-ference"
  end

  # NOTE copy of data from new tree in to replace outisde oldtree
  oldtree.bt = tree.bt
  oldtree.btid = tree.btid
  oldtree.cliques = tree.cliques
  oldtree.frontals = tree.frontals
  oldtree.eliminationOrder = tree.eliminationOrder
  oldtree.buildTime = tree.buildTime

  if opt.drawtree && opt.async
    @warn "due to async=true, only keeping task pointer, not stopping the drawtreerate task!  Consider not using .async together with .drawtreerate != 0"
    push!(smtasks, treetask)
  else
    dotreedraw[1] = 0
  end

  # if debugging and not async then also print the CSMHistory
  if opt.dbg && !opt.async
    hists = !opt.async ? fetchCliqHistoryAll!(smtasks) : hist
    printCSMHistorySequential(hists, joinLogPath(dfgl, "HistoryCSMAll.txt"))
  end

  return oldtree
end

"""
    solveGrapn!

Just an alias, see documentation for `solveTree!`.
"""
DFG.solveGraph!(dfg::AbstractDFG, w...;kw...) = solveTree!(dfg, w...;kw...)

"""
    $SIGNATURES
Internal function used for solveCliqUp! to build the incoming upward message (Rx)
"""
function _buildMessagesUp(
  fg::AbstractDFG,
  tree::AbstractBayesTree,
  cliqid,
  solveKey::Symbol;
  status = UPSOLVED,
)
  #
  cliq = getClique(tree, cliqid)
  beliefMessages = Dict{Int, LikelihoodMessage}()
  for child in getChildren(tree, cliq)
    msg = prepCliqueMsgUp(fg, child, solveKey, status)
    push!(beliefMessages, child.id[] => msg)
  end
  return beliefMessages
end

"""
    $SIGNATURES

Perform inference in the upward direction over one clique in the Bayes tree according to `opt::SolverParams`.

Example
```julia
tree = buildTreeReset!(fg)
hist, upMessageOut = solveCliqUp!(fg, tree, 2)
```

Notes
- Modifies fg with new values
- Calculates up messages from fg if not provided

DevNotes
- Test isfixedlag
- Test recordcliq

Related
[`solveTree!`](@ref), [`buildTreeReset!`](@ref), [`printCliqHistorySummary`](@ref), [`repeatCSMStep!`](@ref), `sandboxStateMachineStep`
"""
function solveCliqUp!(
  fg::AbstractDFG,
  tree::AbstractBayesTree,
  cliqid::Union{CliqueId, Int, Symbol},
  solveKey::Symbol = :default,
  beliefMessages::Dict{Int, LikelihoodMessage} = _buildMessagesUp(
    fg,
    tree,
    cliqid,
    solveKey,
  ); # create belief message from fg if needed
  verbose::Bool = false,
  recordcliq::Bool = false,
)
  # cliqHistories = Dict{Int,Vector{CSMHistoryTuple}}(),
  #

  # hist = Vector{CSMHistoryTuple}()
  opt = DFG.getSolverParams(fg)

  olddown = opt.downsolve
  opt.downsolve = false
  #TODO test 
  if opt.isfixedlag
    @info "Quasi fixed-lag is enabled (a feature currently in testing)!"
    fifoFreeze!(fg)
  end

  cliq = getClique(tree, cliqid)

  # TODO improve, perhaps add to constructor, sommer add all channels here regardless.
  initTreeMessageChannels!(tree)

  @debug "putting messages on up channels from $(keys(beliefMessages))"
  # put the up messages (beliefMessages) that will be used to solve this clique on the channel, the input 
  for (id, msg) in pairs(beliefMessages)
    child = getClique(tree, id)
    for e in getEdgesParent(tree, child)
      @async putBeliefMessageUp!(tree, e, msg)
    end
  end

  #
  @debug "taking belief message that will be sent up"
  # take! the message that is sent up by this clique, the output  
  takeUpTask = @async takeBeliefMessageUp!(tree, getEdgesParent(tree, cliq)[1])

  recordcliqs = recordcliq ? [getFrontals(cliq)[1]] : Symbol[]

  hist = tryCliqStateMachineSolve!(
    fg,
    tree,
    cliq.id;
    solveKey = solveKey,
    verbose = verbose,
    drawtree = opt.drawtree,
    limititers = opt.limititers,
    downsolve = false,
    recordcliqs = recordcliqs,
    incremental = opt.incremental,
  )
  #

  # post-hoc store possible state machine history in clique (without recursively saving earlier history inside state history)
  # assignTreeHistory!(tree, cliqHistories)
  beliefMessageOut = fetch(takeUpTask)
  #restore downsolve
  opt.downsolve = olddown

  return hist, beliefMessageOut
end

"""
    $SIGNATURES
Internal function used for solveCliqDown! to build the incoming downward message (Rx)
"""
function _buildMessageDown(
  fg::AbstractDFG,
  tree::AbstractBayesTree,
  cliqid,
  solveKey::Symbol;
  status::CliqStatus = DOWNSOLVED,
)
  #
  cliq = getClique(tree, cliqid)
  parent = getParent(tree, cliq)[1]
  return getCliqDownMsgsAfterDownSolve(fg, parent, solveKey; status = status)
end

function solveCliqDown!(
  fg::AbstractDFG,
  tree::AbstractBayesTree,
  cliqid::Union{CliqueId, Int, Symbol},
  solveKey::Symbol = :default,
  beliefMessage::LikelihoodMessage = _buildMessageDown(fg, tree, cliqid, solveKey); # create belief message from fg if needed
  verbose::Bool = false,
  recordcliq::Bool = false,
)
  #

  # hist = Vector{CSMHistoryTuple}()
  opt = DFG.getSolverParams(fg)

  upsolve = opt.upsolve

  opt.upsolve = false

  cliq = getClique(tree, cliqid)

  # TODO improve, perhaps add to constructor, sommer add all channels here regardless.
  initTreeMessageChannels!(tree)

  # Build the cliq up message to populate message factors that is needed for down
  @debug "Putting message on up channel from children"
  for (id, msg) in _buildMessagesUp(fg, tree, cliqid, solveKey)
    child = getClique(tree, id)
    for e in getEdgesParent(tree, child)
      @async putBeliefMessageUp!(tree, e, msg)
    end
  end

  # put the down message (beliefMessage) that will be used to solve this clique on the channel, the input 
  @debug "putting message on down channel from parent, used by this clique"
  for e in getEdgesParent(tree, cliq)
    @async putBeliefMessageDown!(tree, e, beliefMessage)
  end

  #take! and discart the up message sent in the skip up part of the solve
  @debug "taking belief message that will be sent up"
  @async takeBeliefMessageUp!(tree, getEdgesParent(tree, cliq)[1])

  #
  @debug "taking belief message that will be sent down"
  # take! the message that is sent down by this clique, the output 
  takeDownTask = @async begin
    messages = Dict{Int, LikelihoodMessage}()
    for e in getEdgesChildren(tree, cliq)
      messages[e.dst] = takeBeliefMessageDown!(tree, e)
    end
    messages
  end

  recordcliqs = recordcliq ? [getFrontals(cliq)[1]] : Symbol[]

  hist = tryCliqStateMachineSolve!(
    fg,
    tree,
    cliq.id;
    solveKey = solveKey,
    verbose = verbose,
    drawtree = opt.drawtree,
    limititers = opt.limititers,
    recordcliqs = recordcliqs,
    incremental = opt.incremental,
  )

  # fetch on down                                  
  beliefMessageOut = fetch(takeDownTask)

  #restore 
  opt.upsolve = upsolve

  return hist, beliefMessageOut
end
