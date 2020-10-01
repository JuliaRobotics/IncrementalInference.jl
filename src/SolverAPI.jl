## Various solver API's used in the past.  These functions are due to be standardized, and obsolete code / functions removed.

export solveTree!, solveGraph!
export fetchCliqHistoryAll!




## ==============================================================================================
# Launch the FSM
## ==============================================================================================

"""
    $SIGNATURES

Perform upward inference using a state machine solution approach.

Notes:
- will call on values from children or parent cliques
- can be called multiple times
- Assumes all cliques in tree are being solved simultaneously and in similar manner.
- State machine rev.1 -- copied from first TreeBasedInitialization.jl.
- Doesn't do partial initialized state properly yet.
"""
function cliqInitSolveUpByStateMachine!(dfg::G,
                                        tree::AbstractBayesTree,
                                        cliq::TreeClique,
                                        timeout::Union{Nothing, <:Real}=nothing;
                                        N::Int=100,
                                        verbose::Bool=false,
                                        verbosefid=stdout,
                                        oldcliqdata::BayesTreeNodeData=BayesTreeNodeData(),
                                        drawtree::Bool=false,
                                        show::Bool=false,
                                        incremental::Bool=true,
                                        limititers::Int=-1,
                                        upsolve::Bool=true,
                                        downsolve::Bool=true,
                                        recordhistory::Bool=false,
                                        delay::Bool=false,
                                        injectDelayBefore::Union{Nothing,Pair{<:Function, <:Real}}=nothing,
                                        logger::SimpleLogger=SimpleLogger(Base.stdout)) where {G <: AbstractDFG, AL <: AbstractLogger}
  #
  children = getChildren(tree, cliq)#Graphs.out_neighbors(cliq, tree.bt)

  prnt = getParent(tree, cliq)

  destType = (G <: InMemoryDFGTypes) ? G : InMemDFGType

  csmc = CliqStateMachineContainer(dfg, initfg(destType, solverParams=getSolverParams(dfg)), tree, cliq, prnt, children, incremental, drawtree, downsolve, delay, getSolverParams(dfg), Dict{Symbol,String}(), oldcliqdata, logger)

  nxt = upsolve ? canCliqMargRecycle_StateMachine : (downsolve ? canCliqMargRecycle_StateMachine : error("must attempt either up or down solve"))

  csmiter_cb = getSolverParams(dfg).drawCSMIters ? ((st::StateMachine)->(cliq.attributes["xlabel"] = st.iter)) : ((st)->())

  statemachine = StateMachine{CliqStateMachineContainer}(next=nxt, name="cliq$(cliq.index)")

  # store statemachine and csmc in task
  if dfg.solverParams.dbg || recordhistory
    task_local_storage(:statemachine, statemachine)
    task_local_storage(:csmc, csmc)
  end

  while statemachine(csmc, timeout, verbose=verbose, verbosefid=verbosefid, verboseXtra=getCliqueStatus(csmc.cliq), iterlimit=limititers, recordhistory=recordhistory, housekeeping_cb=csmiter_cb, injectDelayBefore=injectDelayBefore); end
  statemachine.history
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




function tryCliqStateMachineSolve!( dfg::AbstractDFG,
                                    treel::AbstractBayesTree,
                                    i::Int,
                                    timeout::Union{Nothing, <:Real}=nothing;
                                    verbose::Bool=false,
                                    verbosefid=stdout,
                                    N::Int=100,
                                    oldtree::AbstractBayesTree=emptyBayesTree(),
                                    drawtree::Bool=false,
                                    limititers::Int=-1,
                                    downsolve::Bool=false,
                                    incremental::Bool=false,
                                    injectDelayBefore::Union{Nothing,<:Pair{<:Function, <:Real}}=nothing,
                                    delaycliqs::Vector{Symbol}=Symbol[],
                                    recordcliqs::Vector{Symbol}=Symbol[])
  #
  clst = :na
  cliq = getClique(treel, i)
  syms = getCliqFrontalVarIds(cliq) # ids =
  oldcliq = attemptTreeSimilarClique(oldtree, getCliqueData(cliq))
  oldcliqdata = getCliqueData(oldcliq)
  opts = getSolverParams(dfg)
  # Base.rm(joinpath(opts.logpath,"logs/cliq$i"), recursive=true, force=true)
  mkpath(joinpath(opts.logpath,"logs/cliq$i/"))
  logger = SimpleLogger(open(joinpath(opts.logpath,"logs/cliq$i/log.txt"), "w+")) # NullLogger()
  # global_logger(logger)
  history = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}()
  recordthiscliq = length(intersect(recordcliqs,syms)) > 0
  delaythiscliq = length(intersect(delaycliqs,syms)) > 0
  try
    history = cliqInitSolveUpByStateMachine!(dfg, treel, cliq, timeout, N=N,
                                              verbose=verbose, verbosefid=verbosefid, drawtree=drawtree,
                                              oldcliqdata=oldcliqdata,
                                              injectDelayBefore=injectDelayBefore,
                                              limititers=limititers, downsolve=downsolve, recordhistory=recordthiscliq, incremental=incremental, delay=delaythiscliq, logger=logger )
    #
    if getSolverParams(dfg).dbg || length(history) >= limititers && limititers != -1
      @info "writing logs/cliq$i/csm.txt"
      # @save "/tmp/cliqHistories/cliq$i.jld2" history
      fid = open(joinpath(opts.logpath,"logs/cliq$i/csm.txt"), "w")
      printCliqHistorySummary(fid, history)
      close(fid)
    end
    flush(logger.stream)
    close(logger.stream)
    # clst = getCliqueStatus(cliq)
    # clst = cliqInitSolveUp!(dfg, treel, cliq, drawtree=drawtree, limititers=limititers )
  catch err
    ## TODO -- use this format instead
    # io = IOBuffer()
    # showerror(io, ex, catch_backtrace())
    # err = String(take!(io))
    # msg = "Error while packing '$(f.label)' as '$fnctype', please check the unpacking/packing converters for this factor - \r\n$err"
    # error(msg)

    ## OLD format
    bt = catch_backtrace()
    println()
    showerror(stderr, err, bt)
    # @warn "writing /tmp/caesar/logs/cliq$i/*.txt"
    fid = open(joinpath(opts.logpath,"logs/cliq$i/stacktrace.txt"), "w")
    showerror(fid, err, bt)
    close(fid)
    fid = open(joinpath(opts.logpath,"logs/cliq$(i)_stacktrace.txt"), "w")
    showerror(fid, err, bt)
    close(fid)
    # @save "/tmp/cliqHistories/$(cliq.label).jld2" history
    fid = open(joinpath(opts.logpath,"logs/cliq$i/csm.txt"), "w")
    printCliqHistorySummary(fid, history)
    close(fid)
    fid = open(joinpath(opts.logpath,"logs/cliq$(i)_csm.txt"), "w")
    printCliqHistorySummary(fid, history)
    close(fid)
    flush(logger.stream)
    close(logger.stream)
    error(err)
  end
  # if !(clst in [:upsolved; :downsolved; :marginalized])
  #   error("Clique $(cliq.index), initInferTreeUp! -- cliqInitSolveUp! did not arrive at the desired solution statu: $clst")
  # end
  return history
end


"""
    $SIGNATURES

Fetch solver history from clique state machines that have completed their async Tasks and store in the `hist::Dict{Int,Tuple}` dictionary.
"""
function fetchCliqHistoryAll!(smt::Vector{Task},
                              hist::Dict{Int,Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}=Dict{Int,Vector{Tuple{DateTime, Int,
                                                                    Function, CliqStateMachineContainer}}}() )
  #
  for i in 1:length(smt)
    sm = smt[i]
    # only fetch states that have completed processing
    if sm.state == :done
      haskey(hist, i) ? @warn("overwriting existing history key $i") : nothing
      hist[i] = fetch(sm)
    elseif !isnothing(sm.storage) && haskey(sm.storage, :statemachine)
      hist[i] = sm.storage[:statemachine].history
    end
  end
  hist
end


"""
    $SIGNATURES

Perform tree based initialization of all variables not yet initialized in factor graph as non-blocking method.

Notes:
- To simplify debugging, this method does not include the usual `@ sync` around all the state machine async processes.
- Extract the error stack with a `fetch` on the failed process return by this function.

Related

initInferTreeUp!
"""
function asyncTreeInferUp!( dfg::AbstractDFG,
                            treel::AbstractBayesTree,
                            timeout::Union{Nothing, <:Real}=nothing;
                            oldtree::AbstractBayesTree=emptyBayesTree(),
                            verbose::Bool=false,
                            verbosefid=stdout,
                            drawtree::Bool=false,
                            N::Int=100,
                            limititers::Int=-1,
                            downsolve::Bool=false,
                            incremental::Bool=false,
                            limititercliqs::Vector{Pair{Symbol, Int}}=Pair{Symbol, Int}[],
                            injectDelayBefore::Union{Nothing,Vector{<:Pair{Int,<:Pair{<:Function,<:Real}}}}=nothing,
                            skipcliqids::Vector{Symbol}=Symbol[],
                            delaycliqs::Vector{Symbol}=Symbol[],
                            recordcliqs::Vector{Symbol}=Symbol[] )
  #
  resetTreeCliquesForUpSolve!(treel)
  if drawtree
    pdfpath = joinLogPath(dfg,"bt.pdf")
    drawTree(treel, show=false, filepath=pdfpath)
  end

  # queue all the tasks
  alltasks = Vector{Task}(undef, length(getCliques(treel)))
  # cliqHistories = Dict{Int,Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}()
  if !isTreeSolved(treel, skipinitialized=true)
    # @sync begin
      # duplicate int i into async (important for concurrency)
      for i in 1:length(getCliques(treel))
        scsym = getCliqFrontalVarIds(getClique(treel, i))
        if length(intersect(scsym, skipcliqids)) == 0
          limthiscsm = filter(x -> (x[1] in scsym), limititercliqs)
          limiter = 0<length(limthiscsm) ? limthiscsm[1][2] : limititers
          injDelay = if injectDelayBefore === nothing
            nothing
          else 
            idb = filter((x)->x[1]==i,injectDelayBefore)
            length(idb) == 1 ? idb[1][2] : nothing
          end
          alltasks[i] = @async tryCliqStateMachineSolve!(dfg, treel, i, timeout, oldtree=oldtree, verbose=verbose, verbosefid=verbosefid, drawtree=drawtree, limititers=limiter, downsolve=downsolve, delaycliqs=delaycliqs, recordcliqs=recordcliqs, injectDelayBefore=injDelay, incremental=incremental, N=N)
        end # if
      end # for
    # end # sync
  end # if

  return alltasks #, cliqHistories
end



"""
    $SIGNATURES

Perform tree based initialization of all variables not yet initialized in factor graph.

Related

asyncTreeInferUp!
"""
function initInferTreeUp!(dfg::AbstractDFG,
                          treel::AbstractBayesTree,
                          timeout::Union{Nothing, <:Real}=nothing;
                          oldtree::AbstractBayesTree=emptyBayesTree(),
                          verbose::Bool=false,
                          verbosefid=stdout,
                          drawtree::Bool=false,
                          N::Int=100,
                          limititers::Int=-1,
                          downsolve::Bool=false,
                          incremental::Bool=false,
                          limititercliqs::Vector{Pair{Symbol, Int}}=Pair{Symbol, Int}[],
                          injectDelayBefore::Union{Nothing,Vector{<:Pair{Int,<:Pair{<:Function,<:Real}}}}=nothing,
                          skipcliqids::Vector{Symbol}=Symbol[],
                          recordcliqs::Vector{Symbol}=Symbol[],
                          delaycliqs::Vector{Symbol}=Symbol[],
                          alltasks::Vector{Task}=Task[],
                          runtaskmonitor::Bool=true)
  #
  # revert :downsolved status to :initialized in preparation for new upsolve
  resetTreeCliquesForUpSolve!(treel)
  if drawtree
    pdfpath = joinLogPath(dfg,"bt.pdf")
    drawTree(treel, show=false, filepath=pdfpath)
  end

  # queue all the tasks
  resize!(alltasks,length(getCliques(treel)))
  cliqHistories = Dict{Int,Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}()
  if !isTreeSolved(treel, skipinitialized=true)
    @sync begin
      runtaskmonitor ? (global monitortask = monitorCSMs(treel, alltasks; forceIntExc = true)) : nothing
      # duplicate int i into async (important for concurrency)
      for i in 1:length(getCliques(treel))
        scsym = getCliqFrontalVarIds(getClique(treel, i))
        if length(intersect(scsym, skipcliqids)) == 0
          limthiscsm = filter(x -> (x[1] in scsym), limititercliqs)
          limiter = 0<length(limthiscsm) ? limthiscsm[1][2] : limititers
          injDelay = if injectDelayBefore === nothing
            nothing
          else 
            idb = filter((x)->x[1]==i,injectDelayBefore)
            length(idb) == 1 ? idb[1][2] : nothing
          end
          alltasks[i] = @async tryCliqStateMachineSolve!(dfg, treel, i, timeout, oldtree=oldtree, verbose=verbose, verbosefid=verbosefid, drawtree=drawtree, limititers=limiter, downsolve=downsolve, incremental=incremental, delaycliqs=delaycliqs, injectDelayBefore=injDelay, recordcliqs=recordcliqs,  N=N)
        end # if
      end # for
    end # sync
  end # if

  # if record cliques is in use, else skip computational delay
  0 == length(recordcliqs) ? nothing : fetchCliqHistoryAll!(alltasks, cliqHistories)

  return alltasks, cliqHistories
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
                    multithread=false)
  #
  # workaround in case isolated variables occur
  ensureSolvable!(dfgl)
  opt = getSolverParams(dfgl)

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
  
  hist = Dict{Int, Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}()

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

  # if desired, drawtree in a loop
  treetask, _dotreedraw = drawTreeAsyncLoop(tree, opt; dotreedraw = dotreedraw)

  @info "Do tree based init-inference on tree"
  
  # choose algorithm 
  if algorithm == :parametric
    @error "Under development, do not use, see #539"
    storeOld && @error("parametric storeOld keyword not wired up yet.") 
    initTreeMessageChannels!(tree)
    alltasks, hist = taskSolveTreeParametric!(dfgl, tree; smtasks=smtasks, oldtree=tree, verbose=verbose, drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs, multithread=multithread )
    @info "Finished tree based Parametric inference"
  else #fall back is :default 
    if opt.async
      smtasks = asyncTreeInferUp!(dfgl, tree, timeout, oldtree=oldtree, N=opt.N, verbose=verbose, verbosefid=verbosefid, drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs, limititercliqs=limititercliqs, injectDelayBefore=injectDelayBefore )
      @info "Tree based init-inference progressing asynchronously, check all CSM clique tasks for completion."
    else

      smtasks, hist = initInferTreeUp!(dfgl, tree, timeout; alltasks=smtasks, oldtree=oldtree, N=opt.N, verbose=verbose, verbosefid=verbosefid, drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs, limititercliqs=limititercliqs, injectDelayBefore=injectDelayBefore, runtaskmonitor=runtaskmonitor)

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
    printCliqHistorySequential(hist, nothing, joinLogPath(dfgl,"HistoryCSMAll.txt") )
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
                    # cliqHistories = Dict{Int,Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}(),
                    async::Bool=false )
  #
  # hist = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}()
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
