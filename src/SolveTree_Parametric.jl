"""
    $SIGNATURES

Start tasks (@async or Threads.@spawn threads if multithread=true) to solve (parametric) the factor graph on the tree.
"""
function taskSolveTreeParametric!(dfg::AbstractDFG,
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
                          smtasks=Task[])
  #
  # revert :downsolved status to :initialized in preparation for new upsolve
  #TODO JT die ene lyk reg
  resetTreeCliquesForUpSolve!(treel)
  #TODO JT needs to be updated
  # setTreeCliquesMarginalized!(dfg, treel)

  drawtree ? drawTree(treel, show=true, filepath=joinpath(getSolverParams(dfg).logpath,"bt.pdf")) : nothing

  # queue all the tasks/threads
  resize!(smtasks, getNumCliqs(treel))
  cliqHistories = Dict{Int,Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}()
  if !isTreeSolved(treel, skipinitialized=true)
    @sync begin
      monitortask = monitorCSMs(treel, smtasks)
      # duplicate int i into async (important for concurrency)
      for i in 1:getNumCliqs(treel) # TODO, this might not always work for Graphs.jl
        scsym = getCliqFrontalVarIds(getClique(treel, i))
        if length(intersect(scsym, skipcliqids)) == 0
          if multithread
            smtasks[i] = Threads.@spawn tryCliqStateMachineSolveParametric!(dfg, treel, i, oldtree=oldtree, verbose=verbose, drawtree=drawtree, limititers=limititers, downsolve=downsolve, incremental=incremental, delaycliqs=delaycliqs, recordcliqs=recordcliqs)
          else
            smtasks[i] = @async tryCliqStateMachineSolveParametric!(dfg, treel, i, oldtree=oldtree, verbose=verbose, drawtree=drawtree, limititers=limititers, downsolve=downsolve, incremental=incremental, delaycliqs=delaycliqs, recordcliqs=recordcliqs)
          end
        end # if
      end # for
    end # sync
  end # if

  # if record cliques is in use, else skip computational delay
  0 == length(recordcliqs) ? nothing : fetchCliqHistoryAll!(smtasks, cliqHistories)

  return smtasks, cliqHistories
end

function monitorCSMs(tree, alltasks; forceIntExc::Bool=false)
  task = @async begin
    while true
      all(istaskdone.(alltasks)) && (@info "monitorCSMs: all tasks done"; break)
      for (i,t) in enumerate(alltasks)
        if istaskfailed(t)
          if forceIntExc
            @error "Task $i failed, sending InterruptExceptions to all running CSM tasks"
            throwIntExcToAllTasks(alltasks)
            @debug "done with throwIntExcToAllTasks"
          else
            @error "Task $i failed, sending error to all cliques"
            bruteForcePushErrorCSM(tree)
            # for tree.messages
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
  for (i,t) in enumerate(alltasks)
    if !istaskdone(alltasks[i])
      @debug "Sending InterruptExceptions to CSM task $i"
      schedule(alltasks[i], InterruptException(), error=true)
      @debug "InterruptExceptions CSM task $i"
    end
  end 
  return nothing
end

function bruteForcePushErrorCSM(tree)
    errMsg = LikelihoodMessage(status=:ERROR_STATUS)
    for (i, ch) in tree.messages

        if isready(ch.upMsg)
            take!(ch.upMsg)
        else
            @info("Up edge $i", ch.upMsg)
            @async put!(ch.upMsg, errMsg)
        end
        if isready(ch.downMsg)
            take!(ch.downMsg)
        else
            @info("Down edge $i", ch.downMsg)
            @async put!(ch.downMsg, errMsg)
        end

    end

    for (i, ch) in tree.messages

        while isready(ch.upMsg)
            @info "cleanup take on $i up"
            take!(ch.upMsg)
        end
        while isready(ch.downMsg)
            @info "cleanup take on $i down"
            take!(ch.downMsg)
        end

    end

end

function tryCliqStateMachineSolveParametric!(dfg::G,
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
    history = initStartCliqStateMachineParametric!(dfg, treel, cliq, cliqKey,
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
    error(err)
  end
  # if !(clst in [:upsolved; :downsolved; :marginalized])
  #   error("Clique $(cliq.index), initInferTreeUp! -- cliqInitSolveUp! did not arrive at the desired solution statu: $clst")
  # end
  return history
end
