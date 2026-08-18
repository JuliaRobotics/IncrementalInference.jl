## Various solver API's used in the past.  These functions are due to be standardized, and obsolete code / functions removed.

export solveTree!, solveGraph!
export fetchCliqHistoryAll!

## ==============================================================================================
## Launch the tasks/treads for cliques
## ==============================================================================================

"""
    $TYPEDEF

A `ProgressMeter` meter that every clique state machine may update, serialized by an explicit lock to prevent 
garbled progress spam when multiple tasks try to update it concurrently (not needed for multithreaded mode).
"""
struct LockedProgress{P}
  meter::P
  lock::ReentrantLock
end

LockedProgress(meter) = LockedProgress(meter, ReentrantLock())

ProgressMeter.next!(p::LockedProgress; kw...) = lock(() -> next!(p.meter; kw...), p.lock)
ProgressMeter.finish!(p::LockedProgress; kw...) = lock(() -> finish!(p.meter; kw...), p.lock)

"""
    $SIGNATURES

Start tasks (@async or Threads.@spawn threads if multithread=true) to solve the factor graph on the tree.
"""
function taskSolveTree!(
  dfg::AbstractDFG,
  treel::AbstractBayesTree,
  timeout::Union{Nothing, <:Real} = nothing;
  oldtree::AbstractBayesTree = BayesTree(),
  smtasks = Task[],
  csmoptions::CSMOptions = CSMOptions(;
    solverparams = getSolverParams(dfg),
  ),
)
  #
  # revert DOWNSOLVED status to INITIALIZED in preparation for new upsolve
  resetTreeCliquesForUpSolve!(treel)

  csmoptions.drawtree ? drawTree(treel; show = false, filepath = joinLogPath(dfg, "bt.dot")) : nothing

  cliqHistories = Dict{Int, Vector{CSMHistoryTuple}}()

  resize!(smtasks, getNumCliqs(treel))

  # `ProgressMeter` rewrites one line in place on a terminal, but emits a *new* line per update when
  # stdout is not a TTY (piped to a file, a CI log, a captured `@time`).  The parametric CSM loops
  # internally, so that is thousands of lines of noise — only show the meter where it can redraw.
  approx_iters = _approxCSMIters(getNumCliqs(treel), csmoptions.solver)
  csmoptions.solve_progressbar = if csmoptions.verbose || !isa(stdout, Base.TTY)
    nothing
  else
    LockedProgress(
      ProgressUnknown(; dt = 1, desc = "Solve Progress: approx max $approx_iters, at iter"),
    )
  end

  # queue all the tasks/threads
  if !isTreeSolved(treel; skipinitialized = true)
    @sync begin
      monitortask = monitorCSMs(treel, smtasks)
      # duplicate int i into async (important for concurrency)
      for i = 1:getNumCliqs(treel) # TODO, this might not always work?
        scsym = getCliqFrontalVarIds(getClique(treel, i))
        if length(intersect(scsym, csmoptions.skipcliqids)) == 0
          limthiscsm = filter(x -> (x[1] in scsym), csmoptions.limititercliqs)
          csmoptions.limititers = 0 < length(limthiscsm) ? limthiscsm[1][2] : csmoptions.limititers

          args = (
            dfg,
            treel,
            i,
            timeout,
          )

          smtasks[i] = if csmoptions.multithread
            Threads.@spawn solveClique!(
              args...;
              oldtree,
              csmoptions,
            )
          else
            @async solveClique!(
              args...;
              oldtree,
              csmoptions,
            )
          end
        end # if
      end # for
    end # sync
  end # if

  # if record cliques is in use, else skip computational delay
  0 == length(csmoptions.recordcliqs) ? nothing : fetchCliqHistoryAll!(smtasks, cliqHistories)

  !isnothing(csmoptions.solve_progressbar) && finish!(csmoptions.solve_progressbar)

  return smtasks, cliqHistories
end



function solveClique!(
  dfg::AbstractDFG,
  treel::AbstractBayesTree,
  cliqKey::Union{Int, CliqueId},
  timeout::Union{Nothing, <:Real} = nothing;
  oldtree::AbstractBayesTree = BayesTree(),
  csmoptions::CSMOptions = CSMOptions(;
    solverparams = getSolverParams(dfg),
  ),
  logger::Any = begin
    mkpath(joinpath(csmoptions.solverparams.logpath, "logs"))
    SimpleLogger(open(joinpath(csmoptions.solverparams.logpath, "logs/cliq_($cliqKey).log"), "w+"))
  end
)
  #
  cliq = getClique(treel, cliqKey)
  syms = getCliqFrontalVarIds(cliq)

  oldcliq = attemptTreeSimilarClique(oldtree, getCliqueData(cliq))
  oldcliqdata = getCliqueData(oldcliq)

  opts = csmoptions.solverparams
  # Base.rm(joinpath(opts.logpath,"logs/cliq$i"), recursive=true, force=true)
  mkpath(joinpath(opts.logpath, "logs/cliq$(cliq.id)/"))
  # global_logger(logger)
  history = Vector{CSMHistoryTuple}()
  csmoptions.recordhistory = length(intersect(csmoptions.recordcliqs, syms)) > 0
  csmoptions.delay = length(intersect(csmoptions.delaycliqs, syms)) > 0

  try
    history = initStartCliqStateMachine!(
      dfg,
      treel,
      cliq,
      timeout;
      oldcliqdata,
      logger,
      csmoptions,
    )
    #
    # cliqHistories[cliqKey] = history
    if length(history) >= csmoptions.limititers && csmoptions.limititers != -1
      @debug "writing $(joinpath(opts.logpath, "logs/cliq$(cliq.id)/csm.txt"))"
      # @save "/tmp/cliqHistories/cliq$(cliq.id).jld2" history
      fid = open(joinpath(opts.logpath, "logs/cliq$(cliq.id)/csm.txt"), "w")
      printCliqHistorySummary(fid, history)
      close(fid)
    end
    flush(logger.stream)
    close(logger.stream)
  catch err
    bt = catch_backtrace()
    println()
    showerror(stderr, err, bt)
    @debug "writing $(joinpath(opts.logpath, "logs/cliq$(cliq.id)/stacktrace.txt"))"
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

Build a Bayes (Junction) tree for `dfgl` and initialize its message channels, ready for
[`solveTreePass!`](@ref).

Separated from [`solveTree!`](@ref) so that a tree can be built once and reused across several
passes.
"""
function buildSolveTree!(
  dfgl::AbstractDFG;
  eliminationOrder::Union{Nothing, Vector{Symbol}} = nothing,
  eliminationConstraints::Vector{Symbol} = Symbol[],
)
  opt = getSolverParams(dfgl)
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

  initTreeMessageChannels!(tree)

  return tree
end

"""
    $SIGNATURES

Run one full up/down pass of the clique state machines over an **existing** `tree`, returning
`(smtasks, hist)`.

Dev note: A tree is safe to reuse across passes: [`taskSolveTree!`](@ref) reverts `DOWNSOLVED` cliques to
`INITIALIZED` before starting, and each clique empties its own `upRx` buffer before taking new child
messages, so nothing accumulates from the previous pass.

See also [`buildSolveTree!`](@ref), [`solveTree!`](@ref), [`solveTreeParametric!`](@ref).
"""
function solveTreePass!(
  dfgl::AbstractDFG,
  tree::AbstractBayesTree;
  smtasks::Vector{Task} = Task[],
  oldtree::AbstractBayesTree = BayesTree(),
  csmoptions::CSMOptions = CSMOptions(; solverparams = getSolverParams(dfgl)),
  # which solver to run, taken from `csmoptions` unless named here
  solver::Union{Nothing, AbstractTreeSolver} = csmoptions.solver,
  algorithm::Symbol = isnothing(solver) ? csmoptions.algorithm : _algorithmLabel(solver),
  solveKey::Symbol = csmoptions.solveKey,
)
  csmoptions.solver = solver
  csmoptions.algorithm = algorithm
  csmoptions.solveKey = solveKey

  opt = csmoptions.solverparams
  hist = Dict{Int, Vector{CSMHistoryTuple}}()

  # if desired, drawtree in a loop.  NOTE re-arm the flag: a previous pass sets it to 0 to stop its
  # own draw task, and this vector is shared across passes.
  csmoptions.dotreedraw[1] = 1
  treetask, _dotreedraw = drawTreeAsyncLoop(tree, opt; dotreedraw = csmoptions.dotreedraw)

  @info "Do tree based init-ference"

  _runtasks() = taskSolveTree!(
    dfgl,
    tree,
    csmoptions.timeout;
    smtasks,
    oldtree,
    csmoptions,
  )

  if opt.async
    @async smtasks, hist = _runtasks()
  else
    smtasks, hist = _runtasks()
    @info "Finished tree based init-ference"
  end

  if opt.drawtree && opt.async
    @warn "due to async=true, only keeping task pointer, not stopping the drawtreerate task!  Consider not using .async together with .drawtreerate != 0"
    push!(smtasks, treetask)
  else
    csmoptions.dotreedraw[1] = 0
  end

  return smtasks, hist
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

Example
```julia
# pass in old `tree` to enable compute recycling -- see online Documentation for more details
tree = solveGraph!(fg [,tree])
```

Notes
- Aliased with `solveGraph!` (legacy `solveTree!`)
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
- For example debug usage see #443.


Related

`solveGraph!`, [`solveCliqUp!`](@ref), [`solveCliqDown!`](@ref), [`buildTreeReset!`](@ref), [`repeatCSMStep`](@ref), [`printCSMHistoryLogical`](@ref)
"""
function DistributedFactorGraphs.solveGraph!(
  dfgl::AbstractDFG,
  oldtree::AbstractBayesTree = BayesTree();
  # tree options
  eliminationOrder::Union{Nothing, Vector{Symbol}} = nothing,
  eliminationConstraints::Vector{Symbol} = Symbol[],
  smtasks::Vector{Task} = Task[],
  # solve/execution options
  csmoptions = CSMOptions(;
    solverparams = getSolverParams(dfgl),
  ),
)
  #
  # `:parametric` is a tangent-space solve: one tree pass is a single Gauss-Newton step, so it has
  # its own entry point that builds the tree once and iterates passes over it.
  if csmoptions.algorithm === :parametric
    tree, _, _ = solveTreeParametric!(
      dfgl;
      solver = @something(csmoptions.solver, TangentSpaceSolver()),
      solveKey = csmoptions.solveKey,
      eliminationOrder,
      eliminationConstraints,
      smtasks,
      csmoptions,
    )
    return tree
  end

  # workaround in case isolated variables occur
  ensureSolvable!(dfgl)

  # showtree should force drawtree
  if csmoptions.solverparams.showtree && !csmoptions.solverparams.drawtree
    @info("Since .showtree=true, also bumping .drawtree=true")
  else
    nothing
  end
  csmoptions.solverparams.drawtree |= csmoptions.solverparams.showtree

  # depcrecation
  # update worker pool incase there are more or less
  setWorkerPool!()
  if csmoptions.solverparams.multiproc && nprocs() == 1
    @info "Setting `.multiproc=false` since `Distributed.nprocs() == 1`"
    csmoptions.solverparams.multiproc = false
  end
  
  # NOTE `:parametric` returned above and does its own parametric graphinit in `solveTreeParametric!`
  if csmoptions.solverparams.graphinit
    @info "Ensure variables are all initialized (graphinit)"
    initAll!(dfgl, csmoptions.solveKey)
  end
  # construct tree
  @info "Solving over the Bayes (Junction) tree."

  if csmoptions.solverparams.isfixedlag
    @info "Quasi fixed-lag is enabled (a feature currently in testing, and ignoring solveKey)!"
    fifoFreeze!(dfgl)
  end

  # perhaps duplicate current value
  if csmoptions.storeOld || csmoptions.solverparams.dbg
    ss = listStates(dfgl) .|> string
    ss_ = ss[occursin.(r"default_", ss)] .|> x -> x[9:end]
    filter!(x -> occursin(r"^\d+$", x), ss_)  # ss_ = ss_[occursin.(r"^\d$",ss_)]
    allk = parse.(Int, ss_)
    nextk = length(allk) == 0 ? 0 : maximum(allk) + 1
    newKey = Symbol(:default_, nextk)
    # DFG.cloneStates!(dfgl, newKey, :default; whereSolvable = >=(1))
    for vlabel in ls(dfgl; whereSolvable = >=(1))
      if hasState(dfgl, vlabel, :default)
        DFG.copytoState!(dfgl, vlabel, newKey, getState(dfgl, vlabel, :default))
      end
    end
    # foreach(x->updateVariableSolverData!(dfgl, x, getState(getVariable(dfgl,x), :default), newKey, true, Symbol[]), ls(dfgl, solvable=1))
    @info "storeOld=true, previous :default deepcopied into $newKey for solvable==1 variables."
  end

  !csmoptions.storeOld ? nothing : @error("storeOld keyword not wired up yet.")

  tree = buildSolveTree!(dfgl; eliminationOrder, eliminationConstraints)

  smtasks, hist = solveTreePass!(
    dfgl,
    tree;
    smtasks,
    oldtree,
    csmoptions,
  )

  # NOTE copy of data from new tree in to replace outisde oldtree
  oldtree.bt = tree.bt
  oldtree.btid = tree.btid
  oldtree.cliques = tree.cliques
  oldtree.frontals = tree.frontals
  oldtree.eliminationOrder = tree.eliminationOrder
  oldtree.buildTime = tree.buildTime

  # if debugging and not async then also print the CSMHistory
  if csmoptions.solverparams.dbg && !csmoptions.solverparams.async
    hists = !csmoptions.solverparams.async ? fetchCliqHistoryAll!(smtasks) : hist
    printCSMHistorySequential(hists, joinLogPath(dfgl, "HistoryCSMAll.txt"))
  end

  return oldtree
end


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

  csmoptions = CSMOptions(;
    solverparams = getSolverParams(fg),
    solveKey,
    verbose,
    recordcliqs,
    downsolve = false,
    drawtree = opt.drawtree,
    limititers = opt.limititers,
    incremental = opt.incremental,
  )

  hist = solveClique!(
    fg,
    tree,
    cliq.id;
    csmoptions,
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

  csmoptions = CSMOptions(;
    solverparams = getSolverParams(fg),
    solveKey,
    verbose,
    recordcliqs,
    drawtree = opt.drawtree,
    limititers = opt.limititers,
    incremental = opt.incremental,
  )

  hist = solveClique!(
    fg,
    tree,
    cliq.id;
    csmoptions,
  )

  # fetch on down                                  
  beliefMessageOut = fetch(takeDownTask)

  #restore 
  opt.upsolve = upsolve

  return hist, beliefMessageOut
end
