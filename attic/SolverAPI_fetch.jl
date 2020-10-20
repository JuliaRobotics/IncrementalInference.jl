
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
# error(err)
rethrow()
end
# if !(clst in [:upsolved; :downsolved; :marginalized])
#   error("Clique $(cliq.index), initInferTreeUp! -- cliqInitSolveUp! did not arrive at the desired solution statu: $clst")
# end
return history
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