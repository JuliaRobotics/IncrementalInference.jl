
function taskSolveTree!(dfg::AbstractDFG,
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
  resetTreeCliquesForUpSolve!(treel)

  drawtree ? drawTree(treel, show=true, filepath=joinpath(getSolverParams(dfg).logpath,"bt.pdf")) : nothing

  # queue all the tasks/threads

  cliqHistories = Dict{Int,Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}()
  
  resize!(alltasks, getNumCliqs(treel))

  approx_iters = getNumCliqs(treel)*20
  solve_progressbar = ProgressUnknown("Solve Progress: approx max $approx_iters, at iter")

  if !isTreeSolved(treel, skipinitialized=true)
    @sync begin
      monitortask = IIF.monitorCSMs(treel, alltasks)
      # duplicate int i into async (important for concurrency)
      for i in 1:getNumCliqs(treel) # TODO, this might not always work for Graphs.jl
        scsym = getCliqFrontalVarIds(getClique(treel, i))
        if length(intersect(scsym, skipcliqids)) == 0
          if multithread
            alltasks[i] = Threads.@spawn tryCliqStateMachineSolve_X!(dfg, treel, i, oldtree=oldtree, verbose=verbose, drawtree=drawtree, limititers=limititers, downsolve=downsolve, incremental=incremental, delaycliqs=delaycliqs, recordcliqs=recordcliqs, solve_progressbar=solve_progressbar)
          else
            alltasks[i] = @async tryCliqStateMachineSolve_X!(dfg, treel, i, oldtree=oldtree, verbose=verbose, drawtree=drawtree, limititers=limititers, downsolve=downsolve, incremental=incremental, delaycliqs=delaycliqs, recordcliqs=recordcliqs, solve_progressbar=solve_progressbar)
          end
        end # if
      end # for
    end # sync
  end # if

  # if record cliques is in use, else skip computational delay
  0 == length(recordcliqs) ? nothing : fetchCliqHistoryAll!(alltasks, cliqHistories)
  
  finish!(solve_progressbar)

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
                                             recordcliqs::Vector{Symbol}=Symbol[],
                                             solve_progressbar=nothing) where G <: AbstractDFG
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
                                                    delay=delaythiscliq, logger=logger, solve_progressbar=solve_progressbar )
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
