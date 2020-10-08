
function solveTree_X!(dfgl::AbstractDFG,
                    oldtree::AbstractBayesTree=emptyBayesTree();
                    storeOld::Bool=false,
                    verbose::Bool=false,
                    delaycliqs::Vector{Symbol}=Symbol[],
                    recordcliqs::Vector{Symbol}=Symbol[],
                    limititercliqs::Vector{Pair{Symbol, Int}}=Pair{Symbol, Int}[],
                    skipcliqids::Vector{Symbol}=Symbol[],
                    maxparallel::Union{Nothing, Int}=nothing,
                    variableOrder::Union{Nothing, Vector{Symbol}}=nothing,
                    variableConstraints::Vector{Symbol}=Symbol[],  
                    smtasks::Vector{Task}=Task[],
                    runtaskmonitor::Bool=true,
                    multithread::Bool=false)
  #
  global dotreedraw
  # workaround in case isolated variables occur
  ensureSolvable!(dfgl)
  opt = getSolverParams(dfgl)
  
  #XXX solveTree_X only works with true
  opt.useMsgLikelihoods = true

  # depcrecation
  if maxparallel !== nothing
    @warn "maxparallel keyword is deprecated, use getSolverParams(fg).maxincidence instead."
    opt.maxincidence = maxparallel
  end

  # update worker pool incase there are more or less
  setWorkerPool!()
  if opt.multiproc && nprocs() == 1
    @warn "Cannot use multiproc with only one process, setting `.multiproc=false`."
    opt.multiproc = false
  end

  if opt.graphinit
    @info "ensure all initialized (using graphinit)"
    ensureAllInitialized!(dfgl)
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
  tree = wipeBuildNewTree!(dfgl, variableOrder=variableOrder, drawpdf=opt.drawtree, show=opt.showtree,ensureSolvable=false,filepath=joinpath(opt.logpath,"bt.pdf"), variableConstraints=variableConstraints, ordering=orderMethod)
  # setAllSolveFlags!(tree, false)
  
  initTreeMessageChannels!(tree)
  
  # if desired, drawtree in a loop
  treetask, dotreedraw = drawTreeAsyncLoop(tree, opt)

  @info "Do tree based init-inference on tree"
  if opt.async
    error("not implemented")
    smtasks = asyncTreeInferUp!(dfgl, tree, oldtree=oldtree, N=opt.N, verbose=verbose, drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs, limititercliqs=limititercliqs )
    @info "Tree based init-inference progressing asynchronously, check all CSM clique tasks for completion."
  else
    # smtasks, hist = taskSolveTree_X!(dfgl, tree; alltasks=smtasks, oldtree=oldtree, N=opt.N, verbose=verbose,  drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs, limititercliqs=limititercliqs, runtaskmonitor=runtaskmonitor )
    smtasks, hist = taskSolveTree_X!(dfgl, tree; multithread=multithread, alltasks=smtasks, oldtree=oldtree, verbose=verbose,  drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs )
    @info "Finished tree based init-inference"
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
    #FIXME
    # printCliqHistorySequential(hist, joinLogPath(dfgl,"HistoryCSMAll.txt") )
  end

  return oldtree, smtasks, hist
end


function taskSolveTree_X!(dfg::AbstractDFG,
                          treel::AbstractBayesTree;
                          oldtree::AbstractBayesTree=emptyBayesTree(),
                          drawtree::Bool=false,
                          verbose::Bool=false,
                          limititers::Int=-1,
                          downsolve::Bool=false,
                          incremental::Bool=false,
                          multithread::Bool=false,
                          skipcliqids::Vector{Symbol}=Symbol[],
                          recordcliqs::Vector{Symbol}=Symbol[],
                          delaycliqs::Vector{Symbol}=Symbol[],
                          alltasks = Task[])
  #
  # revert :downsolved status to :initialized in preparation for new upsolve
  #TODO JT die ene lyk reg
  resetTreeCliquesForUpSolve!(treel)
  #TODO JT needs to be updated
  # setTreeCliquesMarginalized!(dfg, treel)

  drawtree ? drawTree(treel, show=true, filepath=joinpath(getSolverParams(dfg).logpath,"bt.pdf")) : nothing

  # queue all the tasks/threads

  cliqHistories = Dict{Int,Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}()
  
  resize!(alltasks, getNumCliqs(treel))

  if !isTreeSolved(treel, skipinitialized=true)
    @sync begin
      monitortask = IIF.monitorCSMs(treel, alltasks)
      # duplicate int i into async (important for concurrency)
      for i in 1:getNumCliqs(treel) # TODO, this might not always work for Graphs.jl
        scsym = getCliqFrontalVarIds(getClique(treel, i))
        if length(intersect(scsym, skipcliqids)) == 0
          if multithread
            alltasks[i] = Threads.@spawn tryCliqStateMachineSolve_X!(dfg, treel, i, oldtree=oldtree, verbose=verbose, drawtree=drawtree, limititers=limititers, downsolve=downsolve, incremental=incremental, delaycliqs=delaycliqs, recordcliqs=recordcliqs)
          else
            alltasks[i] = @async tryCliqStateMachineSolve_X!(dfg, treel, i, oldtree=oldtree, verbose=verbose, drawtree=drawtree, limititers=limititers, downsolve=downsolve, incremental=incremental, delaycliqs=delaycliqs, recordcliqs=recordcliqs)
          end
        end # if
      end # for
    end # sync
  end # if

  # if record cliques is in use, else skip computational delay
  0 == length(recordcliqs) ? nothing : fetchCliqHistoryAll!(alltasks, cliqHistories)

  return alltasks, cliqHistories
end


function tryCliqStateMachineSolve_X!(dfg::G,
                                             treel::AbstractBayesTree,
                                             cliqKey::Int;
                                             oldtree::AbstractBayesTree=emptyBayesTree(),
                                             verbose::Bool=false,
                                             drawtree::Bool=false,
                                             limititers::Int=-1,
                                             downsolve::Bool=false,
                                             incremental::Bool=false,
                                             delaycliqs::Vector{Symbol}=Symbol[],
                                             recordcliqs::Vector{Symbol}=Symbol[]) where G <: AbstractDFG
  #
  clst = :na
  cliq = getClique(treel, cliqKey) #treel.cliques[cliqKey]
  syms = getCliqFrontalVarIds(cliq) # ids =
  # TODO JT Removed old tree reuse
  # oldcliq = attemptTreeSimilarClique(oldtree, getCliqueData(cliq))
  # oldcliqdata = getCliqueData(oldcliq)
  opts = getSolverParams(dfg)
  # Base.rm(joinpath(opts.logpath,"logs/cliq$i"), recursive=true, force=true)
  mkpath(joinpath(opts.logpath,"logs/cliq$(cliq.index)/"))
  logger = SimpleLogger(open(joinpath(opts.logpath,"logs/cliq$(cliq.index)/log.txt"), "w+")) # NullLogger()
  # global_logger(logger)
  history = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}()
  recordthiscliq = length(intersect(recordcliqs,syms)) > 0
  delaythiscliq = length(intersect(delaycliqs,syms)) > 0
  try
    history = initStartCliqStateMachine_X!(dfg, treel, cliq, cliqKey,
                                                    drawtree=drawtree, verbose=verbose,
                                                    limititers=limititers, downsolve=downsolve,
                                                    recordhistory=recordthiscliq, incremental=incremental,
                                                    delay=delaythiscliq, logger=logger )
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
  # if !(clst in [:upsolved; :downsolved; :marginalized])
  #   error("Clique $(cliq.index), initInferTreeUp! -- cliqInitSolveUp! did not arrive at the desired solution statu: $clst")
  # end
  return history
end
