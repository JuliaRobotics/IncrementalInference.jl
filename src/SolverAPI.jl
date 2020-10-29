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
function taskSolveTree!(dfg::AbstractDFG,
                          treel::AbstractBayesTree,
                          timeout::Union{Nothing, <:Real}=nothing;
                          oldtree::AbstractBayesTree=emptyBayesTree(),
                          drawtree::Bool=false,
                          verbose::Bool=false,
                          verbosefid=stdout,
                          limititers::Int=-1,
                          downsolve::Bool=false,
                          incremental::Bool=false,
                          multithread::Bool=false,
                          skipcliqids::Vector{Symbol}=Symbol[],
                          recordcliqs::Vector{Symbol}=Symbol[],
                          delaycliqs::Vector{Symbol}=Symbol[],
                          smtasks=Task[],
                          algorithm::Symbol=:default)
  #
  # revert DOWNSOLVED status to INITIALIZED in preparation for new upsolve
  resetTreeCliquesForUpSolve!(treel)

  drawtree ? drawTree(treel, show=true, filepath=joinpath(getSolverParams(dfg).logpath,"bt.pdf")) : nothing

  cliqHistories = Dict{Int,Vector{CSMHistoryTuple}}()
  
  resize!(smtasks, getNumCliqs(treel))
  
  approx_iters = getNumCliqs(treel)*24
  solve_progressbar = verbose ? nothing : ProgressUnknown("Solve Progress: approx max $approx_iters, at iter")
  
  # queue all the tasks/threads
  if !isTreeSolved(treel, skipinitialized=true)
    @sync begin
      monitortask = monitorCSMs(treel, smtasks)
      # duplicate int i into async (important for concurrency)
      for i in 1:getNumCliqs(treel) # TODO, this might not always work for Graphs.jl
        scsym = getCliqFrontalVarIds(getClique(treel, i))
        if length(intersect(scsym, skipcliqids)) == 0
          if multithread
            smtasks[i] = Threads.@spawn tryCliqStateMachineSolve!(dfg, treel, i, timeout; algorithm=algorithm, oldtree=oldtree, verbose=verbose, verbosefid=verbosefid, drawtree=drawtree, limititers=limititers, downsolve=downsolve, incremental=incremental, delaycliqs=delaycliqs, recordcliqs=recordcliqs, solve_progressbar=solve_progressbar)
          else
            smtasks[i] = @async tryCliqStateMachineSolve!(dfg, treel, i, timeout; algorithm=algorithm, oldtree=oldtree, verbose=verbose, verbosefid=verbosefid, drawtree=drawtree, limititers=limititers, downsolve=downsolve, incremental=incremental, delaycliqs=delaycliqs, recordcliqs=recordcliqs, solve_progressbar=solve_progressbar)
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


function tryCliqStateMachineSolve!(dfg::G,
                                   treel::AbstractBayesTree,
                                   cliqKey::Int,
                                   timeout::Union{Nothing, <:Real}=nothing;
                                   oldtree::AbstractBayesTree=emptyBayesTree(),
                                   verbose::Bool=false,
                                   verbosefid=stdout,
                                   drawtree::Bool=false,
                                   limititers::Int=-1,
                                   downsolve::Bool=false,
                                   incremental::Bool=false,
                                   delaycliqs::Vector{Symbol}=Symbol[],
                                   recordcliqs::Vector{Symbol}=Symbol[],
                                   solve_progressbar=nothing,
                                   algorithm::Symbol=:default) where G <: AbstractDFG
  #
  clst = :na
  cliq = getClique(treel, cliqKey) #treel.cliques[cliqKey]
  syms = getCliqFrontalVarIds(cliq) # ids =
  
  oldcliq = attemptTreeSimilarClique(oldtree, getCliqueData(cliq))
  oldcliqdata = getCliqueData(oldcliq)

  opts = getSolverParams(dfg)
  # Base.rm(joinpath(opts.logpath,"logs/cliq$i"), recursive=true, force=true)
  mkpath(joinpath(opts.logpath,"logs/cliq$(cliq.index)/"))
  logger = SimpleLogger(open(joinpath(opts.logpath,"logs/cliq$(cliq.index)/log.txt"), "w+")) # NullLogger()
  # global_logger(logger)
  history = Vector{CSMHistoryTuple}()
  recordthiscliq = length(intersect(recordcliqs,syms)) > 0
  delaythiscliq = length(intersect(delaycliqs,syms)) > 0
  try
    history = initStartCliqStateMachine!(dfg, treel, cliq, timeout;
                                         oldcliqdata=oldcliqdata,
                                         drawtree=drawtree, verbose=verbose, verbosefid=verbosefid,
                                         limititers=limititers, downsolve=downsolve,
                                         recordhistory=recordthiscliq, incremental=incremental,
                                         delay=delaythiscliq, logger=logger, solve_progressbar=solve_progressbar,
                                         algorithm=algorithm )
    #
    # cliqHistories[cliqKey] = history
    if length(history) >= limititers && limititers != -1
      # @warn "writing logs/cliq$(cliq.index)/csm.txt"
      # @save "/tmp/cliqHistories/cliq$(cliq.index).jld2" history
      fid = open(joinpath(opts.logpath,"logs/cliq$(cliq.index)/csm.txt"), "w")
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
    # @warn "writing /tmp/caesar/logs/cliq$(cliq.index)/*.txt"
    fid = open(joinpath(opts.logpath,"logs/cliq$(cliq.index)/stacktrace.txt"), "w")
    showerror(fid, err, bt)
    close(fid)
    fid = open(joinpath(opts.logpath,"logs/cliq$(cliq.index)_stacktrace.txt"), "w")
    showerror(fid, err, bt)
    close(fid)
    # @save "/tmp/cliqHistories/$(cliq.label).jld2" history
    fid = open(joinpath(opts.logpath,"logs/cliq$(cliq.index)/csm.txt"), "w")
    printCliqHistorySummary(fid, history)
    close(fid)
    fid = open(joinpath(opts.logpath,"logs/cliq$(cliq.index)_csm.txt"), "w")
    printCliqHistorySummary(fid, history)
    close(fid)
    flush(logger.stream)
    close(logger.stream)
    rethrow()
  end
  # if !(clst in [UPSOLVED; DOWNSOLVED; MARGINALIZED])
  #   error("Clique $(cliq.index), initInferTreeUp! -- cliqInitSolveUp! did not arrive at the desired solution statu: $clst")
  # end
  return history
end


"""
    $SIGNATURES

Standalone state machine solution for a single clique.

Related:

initInferTreeUp!
"""
function solveCliqWithStateMachine!(dfg::G,
                                    tree::AbstractBayesTree,
                                    frontal::Symbol;
                                    iters::Int=200,
                                    downsolve::Bool=true,
                                    recordhistory::Bool=false,
                                    verbose::Bool=false,
                                    nextfnc::Function=canCliqMargRecycle_StateMachine,
                                    prevcsmc::Union{Nothing,CliqStateMachineContainer}=nothing) where G <: AbstractDFG
  #
  cliq = getClique(tree, frontal)

  children = getChildren(tree, cliq)#Graphs.out_neighbors(cliq, tree.bt)

  prnt = getParent(tree, cliq)

  destType = (G <: InMemoryDFGTypes) ? G : InMemDFGType

  csmc = isa(prevcsmc, Nothing) ? CliqStateMachineContainer(dfg, initfg(destType, solverParams=getSolverParams(dfg)), tree, cliq, prnt, children, false, true, true, downsolve, false, getSolverParams(dfg)) : prevcsmc
  statemachine = StateMachine{CliqStateMachineContainer}(next=nextfnc, name="cliq$(cliq.index)")
  while statemachine(csmc, verbose=verbose, iterlimit=iters, recordhistory=recordhistory); end
  statemachine, csmc
end


## ==============================================================================================
# Prepare CSM (based on FSM) entry points
## ==============================================================================================



"""
    $SIGNATURES

Fetch solver history from clique state machines that have completed their async Tasks and store in the `hist::Dict{Int,Tuple}` dictionary.
"""
function fetchCliqHistoryAll!(smt::Vector{Task},
                              hist::Dict{Int,Vector{CSMHistoryTuple}}=Dict{Int,Vector{CSMHistoryTuple}}() )
  #
  for i in 1:length(smt)
    sm = smt[i]
    # only fetch states that have completed processing
    if sm.state == :done
      haskey(hist, i) ? @warn("overwriting existing history key $i") : nothing
      hist[i] = fetch(sm)
    elseif !isnothing(sm.storage) && haskey(sm.storage, :statemachine)
      hist[i] = CSMHistoryTuple.(sm.storage[:statemachine].history)
    end
  end
  hist
end






## ==============================================================================================
# Nominal user interface to the solver
## ==============================================================================================




"""
    $SIGNATURES

Perform inference over the Bayes tree according to `opt::SolverParams`.

Notes
- Variety of options, including fixed-lag solving -- see `getSolverParams(fg)` for details.
- Latest result always stored in `solvekey=:default`.
- Experimental `storeOld::Bool=true` will duplicate the current result as supersolve `:default_k`.
  - Based on `solvable==1` assumption.
- `limititercliqs` allows user to limit the number of iterations a specific CSM does.
- keywords `verbose` and `verbosefid::IOStream` can be used together to to send output to file or default `stdout`.

Example
```julia
# without [or with] compute recycling
tree, smt, hist = solveTree!(fg [,tree])
```

Related

solveCliq!, resetBuildTree!
"""
function solveTree!(dfgl::AbstractDFG,
                    oldtree::AbstractBayesTree=emptyBayesTree();
                    timeout::Union{Nothing, <:Real}=nothing,
                    storeOld::Bool=false,
                    verbose::Bool=false,
                    verbosefid=stdout,
                    delaycliqs::Vector{Symbol}=Symbol[],
                    recordcliqs::Vector{Symbol}=Symbol[],
                    limititercliqs::Vector{Pair{Symbol, Int}}=Pair{Symbol, Int}[],
                    injectDelayBefore::Union{Nothing,Vector{<:Pair{Int,<:Pair{<:Function, <:Real}}}}=nothing,
                    skipcliqids::Vector{Symbol}=Symbol[],
                    maxparallel::Union{Nothing, Int}=nothing,
                    variableOrder::Union{Nothing, Vector{Symbol}}=nothing,
                    variableConstraints::Vector{Symbol}=Symbol[],
                    smtasks::Vector{Task}=Task[],
                    dotreedraw = Int[1;],
                    runtaskmonitor::Bool=true,
                    algorithm::Symbol=:default,
                    multithread::Bool=false)
  #
  # workaround in case isolated variables occur
  ensureSolvable!(dfgl)
  opt = getSolverParams(dfgl)

  opt.useMsgLikelihoods == false && @warn("#TODO Verify useMsgLikelihoods=false works correctly")

  # depcrecation
  if maxparallel !== nothing
    @warn "`maxparallel` keyword is deprecated, use `getSolverParams(fg).maxincidence` instead."
    opt.maxincidence = maxparallel
  end

  # update worker pool incase there are more or less
  setWorkerPool!()
  if opt.multiproc && nprocs() == 1
    @info "Setting `.multiproc=false` since `Distributed.nprocs() == 1`"
    opt.multiproc = false
  end

  if opt.graphinit
    @info "Ensure variables are all initialized (graphinit)"
    ensureAllInitialized!(dfgl)
    if algorithm==:parametric
      @warn "Parametric is using default graphinit"
      initParametricFrom(dfgl)
    end
  end
  # construct tree
  @info "Solving over the Bayes (Junction) tree."
  
  hist = Dict{Int, Vector{CSMHistoryTuple}}()

  if opt.isfixedlag
      @info "Quasi fixed-lag is enabled (a feature currently in testing)!"
      fifoFreeze!(dfgl)
  end

  # perhaps duplicate current value
  if storeOld || opt.dbg
    ss = listSupersolves(dfgl) .|> string
    ss_ = ss[occursin.(r"default_",ss)] .|> x->x[9:end]
    filter!(x->occursin(r"^\d+$",x), ss_)  # ss_ = ss_[occursin.(r"^\d$",ss_)]
    allk = parse.(Int,ss_)
    nextk = length(allk) == 0 ? 0 : maximum(allk)+1
    newKey = Symbol(:default_, nextk)
    deepcopySupersolve!(dfgl, newKey, :default, solvable=1)
    # foreach(x->updateVariableSolverData!(dfgl, x, getSolverData(getVariable(dfgl,x), :default), newKey, true, Symbol[]), ls(dfgl, solvable=1))
    @info "storeOld=true, previous :default deepcopied into $newKey for solvable==1 variables."
  end

  orderMethod = 0 < length(variableConstraints) ? :ccolamd : :qr

  # current incremental solver builds a new tree and matches against old tree for recycling.
  tree = resetBuildTree!(dfgl, variableOrder=variableOrder, drawpdf=opt.drawtree, show=opt.showtree,ensureSolvable=false,filepath=joinpath(opt.logpath,"bt.pdf"), variableConstraints=variableConstraints, ordering=orderMethod)
  # setAllSolveFlags!(tree, false)
  
  initTreeMessageChannels!(tree)
  
  # if desired, drawtree in a loop
  treetask, _dotreedraw = drawTreeAsyncLoop(tree, opt; dotreedraw = dotreedraw)

  @info "Do tree based init-inference on tree"

  # choose algorithm 
  if algorithm == :parametric
    @error "Under development, do not use, see #539"
    storeOld && @error("parametric storeOld keyword not wired up yet.") 
    # alltasks, hist = taskSolveTreeParametric!(dfgl, tree; smtasks=smtasks, oldtree=tree, verbose=verbose, drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs, multithread=multithread )
    smtasks, hist = taskSolveTree!(dfgl, tree, timeout; algorithm=algorithm, multithread=multithread, smtasks=smtasks, oldtree=oldtree,          verbose=verbose, verbosefid=verbosefid, drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs)
    @info "Finished tree based Parametric inference"
  else # fall back is :default with take CSM
    if opt.async
      @async smtasks, hist = taskSolveTree!(dfgl, tree, timeout;algorithm=algorithm, multithread=multithread, smtasks=smtasks, oldtree=oldtree,          verbose=verbose,                        drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs)
    else
      # smtasks, hist = taskSolveTree!(dfgl, tree; alltasks=smtasks, oldtree=oldtree, N=opt.N, verbose=verbose,  drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs, limititercliqs=limititercliqs, runtaskmonitor=runtaskmonitor )
      smtasks, hist = taskSolveTree!(dfgl, tree, timeout; algorithm=algorithm, multithread=multithread, smtasks=smtasks, oldtree=oldtree,          verbose=verbose, verbosefid=verbosefid, drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs)
    # smtasks, hist = initInferTreeUp!(dfgl, tree, timeout;                                            alltasks=smtasks, oldtree=oldtree, N=opt.N, verbose=verbose, verbosefid=verbosefid, drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs, limititercliqs=limititercliqs, injectDelayBefore=injectDelayBefore, runtaskmonitor=runtaskmonitor)
      @info "Finished tree based init-inference"
    end
  end


  # NOTE copy of data from new tree in to replace outisde oldtree
  oldtree.bt = tree.bt
  oldtree.btid = tree.btid
  oldtree.cliques = tree.cliques
  oldtree.frontals = tree.frontals
  oldtree.variableOrder = tree.variableOrder
  oldtree.buildTime = tree.buildTime

  hist = !opt.async ? fetchCliqHistoryAll!(smtasks) : hist

  if opt.drawtree && opt.async
    @warn "due to async=true, only keeping task pointer, not stopping the drawtreerate task!  Consider not using .async together with .drawtreerate != 0"
    push!(smtasks, treetask)
  else
    dotreedraw[1] = 0
  end

  # if debugging and not async then also print the CSMHistory
  if opt.dbg && !opt.async
    printCSMHistorySequential(hist, joinLogPath(dfgl,"HistoryCSMAll.txt") )
  end

  return oldtree, smtasks, hist
end

"""
$SIGNATURES
See `solveTree!`.
"""
const solveGraph! = solveTree!

"""
    $SIGNATURES

Perform inference over one clique in the Bayes tree according to `opt::SolverParams`.

Example
```julia
tree = resetBuildTree!(fg)
smt, hist = solveCliq!(fg, tree, :x1 [,cliqHistories=hist] )
```

Related

solveTree!, resetBuildTree!
"""
function solveCliq!(dfgl::AbstractDFG,
                    tree::AbstractBayesTree,
                    cliqid::Symbol;
                    verbose::Bool=false,
                    recordcliq::Bool=false,
                    # cliqHistories = Dict{Int,Vector{CSMHistoryTuple}}(),
                    async::Bool=false )
  #
  # hist = Vector{CSMHistoryTuple}()
  opt = DFG.getSolverParams(dfgl)

  if opt.isfixedlag
      @info "Quasi fixed-lag is enabled (a feature currently in testing)!"
      fifoFreeze!(dfgl)
  end

  # if !isTreeSolved(treel, skipinitialized=true)
  cliq = getClique(tree, cliqid)
  cliqtask = if async
    @async tryCliqStateMachineSolve!(dfgl, tree, cliq.index, verbose=verbose, drawtree=opt.drawtree, limititers=opt.limititers, downsolve=opt.downsolve,recordcliqs=(recordcliq ? [cliqid] : Symbol[]), incremental=opt.incremental)
  else
    tryCliqStateMachineSolve!(dfgl, tree, cliq.index, verbose=verbose, drawtree=opt.drawtree, limititers=opt.limititers, downsolve=opt.downsolve,recordcliqs=(recordcliq ? [cliqid] : Symbol[]), incremental=opt.incremental) # N=N
  end
  # end # if

  # post-hoc store possible state machine history in clique (without recursively saving earlier history inside state history)
  # assignTreeHistory!(tree, cliqHistories)

  # cliqHistories
  return cliqtask
end
