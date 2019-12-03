"""
    $SIGNATURES

Start async tasks to solve (parametric) on the tree.
Related
TODO multithrededSolveTreeParametric!
"""
function taskSolveTreeParametric!(dfg::AbstractDFG,
                          treel::BayesTree;
                          oldtree::BayesTree=emptyBayesTree(),
                          drawtree::Bool=false,
                          limititers::Int=-1,
                          downsolve::Bool=false,
                          incremental::Bool=false,
                          skipcliqids::Vector{Symbol}=Symbol[],
                          recordcliqs::Vector{Symbol}=Symbol[],
                          delaycliqs::Vector{Symbol}=Symbol[])
  #
  # revert :downsolved status to :initialized in preparation for new upsolve
  #TODO JT die ene lyk reg
  resetTreeCliquesForUpSolve!(treel)
  #TODO JT needs to be updated
  # setTreeCliquesMarginalized!(dfg, treel)

  drawtree ? drawTree(treel, show=true, filepath=joinpath(getSolverParams(dfg).logpath,"bt.pdf")) : nothing

  # queue all the tasks
  alltasks = Vector{Task}(undef, length(treel.cliques))
  cliqHistories = Dict{Int,Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}()
  if !isTreeSolved(treel, skipinitialized=true)
    @sync begin
      # duplicate int i into async (important for concurrency)
      for i in 1:length(treel.cliques)
        scsym = getCliqFrontalVarIds(treel.cliques[i])
        if length(intersect(scsym, skipcliqids)) == 0
          alltasks[i] = @async tryCliqStateMachineSolveParametric!(dfg, treel, i, oldtree=oldtree, drawtree=drawtree, limititers=limititers, downsolve=downsolve, incremental=incremental, delaycliqs=delaycliqs, recordcliqs=recordcliqs)
        end # if
      end # for
    end # sync
  end # if

  fetchCliqTaskHistoryAll!(alltasks, cliqHistories)

  # post-hoc store possible state machine history in clique (without recursively saving earlier history inside state history)
  assignTreeHistory!(treel, cliqHistories)
  # for i in 1:length(treel.cliques)
  #   if haskey(cliqHistories, i)
  #     hist = cliqHistories[i]
  #     for i in 1:length(hist)
  #       hist[i][4].logger = ConsoleLogger()
  #     end
  #     getData(treel.cliques[i]).statehistory=hist
  #   end
  # end

  return alltasks, cliqHistories
end


function tryCliqStateMachineSolveParametric!(dfg::G,
                                             treel::BayesTree,
                                             i::Int;
                                             # cliqHistories;
                                             oldtree::BayesTree=emptyBayesTree(),
                                             drawtree::Bool=false,
                                             limititers::Int=-1,
                                             downsolve::Bool=false,
                                             incremental::Bool=false,
                                             delaycliqs::Vector{Symbol}=Symbol[],
                                             recordcliqs::Vector{Symbol}=Symbol[]) where G <: AbstractDFG
  #
  clst = :na
  cliq = treel.cliques[i]
  syms = getCliqFrontalVarIds(cliq) # ids =
  # TODO JT Removed old tree reuse
  # oldcliq = attemptTreeSimilarClique(oldtree, getData(cliq))
  # oldcliqdata = getData(oldcliq)
  opts = getSolverParams(dfg)
  # Base.rm(joinpath(opts.logpath,"logs/cliq$i"), recursive=true, force=true)
  mkpath(joinpath(opts.logpath,"logs/cliq$i/"))
  logger = SimpleLogger(open(joinpath(opts.logpath,"logs/cliq$i/log.txt"), "w+")) # NullLogger()
  # global_logger(logger)
  history = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}()
  recordthiscliq = length(intersect(recordcliqs,syms)) > 0
  delaythiscliq = length(intersect(delaycliqs,syms)) > 0
  try
    history = initStartCliqStateMachineParametric!(dfg, treel, cliq, drawtree=drawtree,
                                             # oldcliqdata=oldcliqdata,
                                             limititers=limititers, downsolve=downsolve, recordhistory=recordthiscliq, incremental=incremental, delay=delaythiscliq, logger=logger )
    #
    # cliqHistories[i] = history
    if length(history) >= limititers && limititers != -1
      # @warn "writing logs/cliq$i/csm.txt"
      # @save "/tmp/cliqHistories/cliq$i.jld2" history
      fid = open(joinpath(opts.logpath,"logs/cliq$i/csm.txt"), "w")
      printCliqHistorySummary(fid, history)
      close(fid)
    end
    flush(logger.stream)
    close(logger.stream)
    # clst = getCliqStatus(cliq)
    # clst = cliqInitSolveUp!(dfg, treel, cliq, drawtree=drawtree, limititers=limititers )
  catch err
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
