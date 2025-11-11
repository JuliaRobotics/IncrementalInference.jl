
# starting to add exports here

global WORKERPOOL = WorkerPool()

"""
    $SIGNATURES

For use with `multiproc`, nominal use is a worker pool of all processes available above and including 2..., but will return single process [1;] if only the first processes is available.
"""
function setWorkerPool!(pool::Vector{Int} = 1 < nprocs() ? setdiff(procs(), [1;]) : [1;])
  global WORKERPOOL
  return WORKERPOOL = WorkerPool(pool)
end

function getWorkerPool()
  global WORKERPOOL
  return WORKERPOOL
end

## =============================================================================================
# Iterate over variables in clique
## =============================================================================================

"""
    $SIGNATURES

Dev Notes
- part of refactoring fmcmc.
- function seems excessive
"""
function compileFMCMessages(
  fgl::AbstractDFG,
  lbls::Vector{Symbol},
  solveKey::Symbol,
  logger = ConsoleLogger(),
)
  #
  d = Dict{Symbol, TreeBelief}()
  for vsym in lbls
    vari = DFG.getVariable(fgl, vsym)
    d[vsym] = TreeBelief(vari, solveKey)
  end
  return d
end

function doFMCIteration(
  fgl::AbstractDFG,
  vsym::Symbol,
  solveKey::Symbol,
  cliq::TreeClique,
  fmsgs,
  N::Int,
  dbg::Bool,
  needFreshMeasurements::Bool = true,
  logger = ConsoleLogger(),
)
  #

  vert = DFG.getVariable(fgl, vsym)
  if !getState(vert, solveKey).marginalized
    # potprod = nothing
    dens, ipc = propagateBelief(
      fgl,
      vsym,
      :;
      needFreshMeasurements,
      N,
      dbg,
      logger,
    )

    if 0 < Npts(dens)
      setBelief!(vert, dens, true, ipc)
      # setValKDE!(vert, densPts, true, ipc)
      # TODO perhaps more debugging inside `propagateBelief`?
    end
  end
  return nothing
end

"""
    $(SIGNATURES)

Iterate successive approximations of clique marginal beliefs by means
of the stipulated proposal convolutions and products of the functional objects
for tree clique `cliq`.
"""
function fmcmc!(
  fgl::AbstractDFG,
  cliq::TreeClique,
  fmsgs::Vector{<:LikelihoodMessage},
  lbls::Vector{Symbol},
  solveKey::Symbol,
  N::Int,
  MCMCIter::Int,
  dbg::Bool = false,
  logger = ConsoleLogger(),
  multithreaded::Bool = false,
)
  #
  with_logger(logger) do
    @info "---------- successive fnc approx ------------$(getLabel(cliq))"
  end
  # repeat several iterations of functional Gibbs sampling for fixed point convergence
  if length(lbls) == 1
    MCMCIter = 1
  end
  # mcmcdbg = Array{CliqGibbsMC,1}()

  # burn-in loop for outer Gibbs
  for iter = 1:MCMCIter
    # iterate through each of the variables, KL-divergence tolerence would be nice test here
    with_logger(logger) do
      @info "#$(iter)\t -- "
    end
    # dbgvals = !dbg ? nothing : CliqGibbsMC([], Symbol[])

    needFreshMeasurements = iter == 1 || getSolverParams(fgl).alwaysFreshMeasurements

    # outer Gibbs cycle
    for vsym in lbls
      doFMCIteration(
        fgl,
        vsym,
        solveKey,
        cliq,
        fmsgs,
        N,
        dbg,
        needFreshMeasurements,
        logger,
      )
    end
    # !dbg ? nothing : push!(mcmcdbg, dbgvals)
  end

  # populate dictionary for return NBPMessage in multiple dispatch
  msgdict = compileFMCMessages(fgl, lbls, solveKey, logger)

  return msgdict
end

## =============================================================================================
# Up solve
## =============================================================================================

"""
    $(SIGNATURES)

Perform computations required for the upward message passing during belief propation on the Bayes (Junction) tree.
This function is usually called as via remote_call for multiprocess dispatch.

Notes
- `fg` factor graph,
- `tree` Bayes tree,
- `cliq` which cliq to perform the computation on,
- `parent` the parent clique to where the upward message will be sent,
- `childmsgs` is for any incoming messages from child cliques.

DevNotes
- FIXME total rewrite with AMP #41 and RoME #244 in mind
"""
function upGibbsCliqueDensity(
  dfg::AbstractDFG,
  cliq::TreeClique,
  solveKey::Symbol,
  inmsgs,
  N::Int = getSolverParams(dfg).N,
  dbg::Bool = false,
  iters::Int = 3,
  logger = ConsoleLogger(),
) # where {T, T2}
  #
  with_logger(logger) do
    @info "up w $(length(inmsgs)) msgs"
  end
  # TODO -- some weirdness with: d,. = d = ., nothing
  # mcmcdbg = Array{CliqGibbsMC,1}()
  d = Dict{Symbol, TreeBelief}()

  # priorprods = Vector{CliqGibbsMC}()

  cliqdata = getCliqueData(cliq)

  # use nested structure for more efficient Chapman-Kolmogorov solution approximation
  if false
    IDS = [cliqdata.frontalIDs; cliqdata.separatorIDs] #inp.cliq.attributes["frontalIDs"]
    d = fmcmc!(dfg, cliq, inmsgs, IDS, solveKey, N, iters, dbg, logger)
  else
    # NOTE -- previous mistake, must iterate over directsvarIDs also (or incorporate once at the right time)
    # NOTE -- double up on directs to allow inflation to take proper affect, see #1051
    d = fmcmc!(
      dfg,
      cliq,
      inmsgs,
      cliqdata.directFrtlMsgIDs,
      solveKey,
      N,
      1,
      dbg,
      logger,
      true,
    )

    if length(cliqdata.msgskipIDs) > 0
      dd = fmcmc!(dfg, cliq, inmsgs, cliqdata.msgskipIDs, solveKey, N, 1, dbg, logger, true)
      for md in dd
        d[md[1]] = md[2]
      end
    end
    if length(cliqdata.itervarIDs) > 0
      ddd = fmcmc!(
        dfg,
        cliq,
        inmsgs,
        cliqdata.itervarIDs,
        solveKey,
        N,
        iters,
        dbg,
        logger,
        false,
      )
      for md in ddd
        d[md[1]] = md[2]
      end
    end
    if length(cliqdata.directPriorMsgIDs) > 0
      doids = setdiff(cliqdata.directPriorMsgIDs, cliqdata.msgskipIDs)
      dddd = fmcmc!(dfg, cliq, inmsgs, doids, solveKey, N, 1, dbg, logger, true)
      for md in dddd
        d[md[1]] = md[2]
      end
    end
  end

  return d
end

## =============================================================================================
# Down solve
## =============================================================================================

## ============================================================================
# Initialization is slightly different and likely to be consolidated
## ============================================================================

"""
    $SIGNATURES

Cycle through var order and initialize variables as possible in `subfg::AbstractDFG`.
Return true if something was updated.

Notes:
- assumed `subfg` is a subgraph containing only the factors that can be used.
  - including the required up or down messages
- intended for both up and down initialization operations.

Dev Notes
- Should monitor updates based on the number of inferred & solvable dimensions
"""
function cycleInitByVarOrder!(
  subfg::AbstractDFG,
  varorder::Vector{Symbol};
  solveKey::Symbol = :default,
  logger = ConsoleLogger(),
)
  #
  with_logger(logger) do
    @info "cycleInitByVarOrder! -- varorder=$(varorder)"
  end
  retval = false
  count = 1
  while count > 0
    count = 0
    for vsym in varorder
      var = DFG.getVariable(subfg, vsym)
      isinit = isInitialized(var, solveKey)
      with_logger(logger) do
        @info "var.label=$(var.label) is initialized=$(isinit)"
      end
      doautoinit!(subfg, [var;]; solveKey = solveKey, logger = logger)
      if isinit != isInitialized(var, solveKey)
        count += 1
        retval = true
      end
    end
  end
  with_logger(logger) do
    @info "cycleInitByVarOrder!, retval=$(retval)"
  end
  flush(logger.stream)
  return retval
end

#
