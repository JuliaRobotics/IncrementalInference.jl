## Various solver API's used in the past.  These functions are due to be standardized, and obsolete code / functions removed.


"""
    $SIGNATURES

Perform inference over the Bayes tree according to `opt::SolverParams`.

Notes
- Variety of options, including fixed-lag solving -- see `getSolverParams(fg)` for details.
- Latest result always stored in `solvekey=:default`.
- Experimental `storeOld::Bool=true` will duplicate the current result as supersolve `:default_k`.
  - Based on `solvable==1` assumption.

Example
```julia
# without [or with] compute recycling
tree, smt, hist = solveTree!(fg [,tree])
```

Related

solveCliq!, wipeBuildNewTree!
"""
function solveTree!(dfgl::G,
                    oldtree::AbstractBayesTree=emptyBayesTree();
                    storeOld::Bool=false,
                    verbose::Bool=false,
                    delaycliqs::Vector{Symbol}=Symbol[],
                    recordcliqs::Vector{Symbol}=Symbol[],
                    skipcliqids::Vector{Symbol}=Symbol[],
                    maxparallel::Union{Nothing, Int}=nothing,
                    variableOrder::Union{Nothing, Vector{Symbol}}=nothing,
                    variableConstraints::Vector{Symbol}=Symbol[]  ) where G <: DFG.AbstractDFG
  #
  # workaround in case isolated variables occur
  ensureSolvable!(dfgl)
  opt = getSolverParams(dfgl)

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
  smtasks=Vector{Task}()
  hist = Dict{Int, Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}()

  if opt.isfixedlag
      @info "Quasi fixed-lag is enabled (a feature currently in testing)!"
      fifoFreeze!(dfgl)
  end

  # perhaps duplicate current value
  if storeOld || opt.dbg
    ss = listSupersolves(dfgl) .|> string
    ss_ = ss[occursin.(r"default_",ss)] .|> x->x[9:end]
    filter!(x->occursin(r"^\d$",x), ss_)  # ss_ = ss_[occursin.(r"^\d$",ss_)]
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

  # if desired, drawtree in a loop
  treetask, dotreedraw = drawTreeAsyncLoop(tree, opt )

  @info "Do tree based init-inference on tree"
  if opt.async
    smtasks = asyncTreeInferUp!(dfgl, tree, oldtree=oldtree, N=opt.N, verbose=verbose, drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs )
    @info "Tree based init-inference progressing asynchronously, check all CSM clique tasks for completion."
  else
    smtasks, hist = initInferTreeUp!(dfgl, tree, oldtree=oldtree, N=opt.N, verbose=verbose,  drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs )
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
    printCliqHistorySequential(hist, joinLogPath(dfgl,"HistoryCSMAll.txt") )
  end

  return oldtree, smtasks, hist
end

"""
$SIGNATURES
See solveTree!.
"""
const solveGraph! = solveTree!

"""
    $SIGNATURES

Perform inference over one clique in the Bayes tree according to `opt::SolverParams`.

Example
```julia
tree = wipeBuildNewTree!(fg)
smt, hist = solveCliq!(fg, tree, :x1 [,cliqHistories=hist] )
```

Related

solveTree!, wipeBuildNewTree!
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



## Experimental Parametric
"""
    $SIGNATURES

Perform parametric inference over the Bayes tree according to `opt::SolverParams`.

Example
```julia
tree, smt, hist = solveTree!(fg ,tree)
```
"""
function solveTreeParametric!(dfgl::DFG.AbstractDFG,
                              tree::AbstractBayesTree;
                              storeOld::Bool=false,
                              verbose::Bool=false,
                              delaycliqs::Vector{Symbol}=Symbol[],
                              recordcliqs::Vector{Symbol}=Symbol[],
                              skipcliqids::Vector{Symbol}=Symbol[],
                              maxparallel::Union{Nothing, Int}=nothing  )
  #
  @error "Under development, do not use, see #539"
  @info "Solving over the Bayes (Junction) tree."
  smtasks=Vector{Task}()
  hist = Dict{Int, Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}()
  opt = DFG.getSolverParams(dfgl)

  # depcrecation
  if maxparallel !== nothing
    @warn "maxparallel keyword is deprecated, use getSolverParams(fg).maxincidence instead."
    opt.maxincidence = maxparallel
  end

  storeOld ? @error("parametric storeOld keyword not wired up yet.") : nothing

  # update worker pool incase there are more or less
  setWorkerPool!()
  if getSolverParams(dfgl).multiproc && nprocs() == 1
    @warn "Cannot use multiproc with only one process, setting `.multiproc=false`."
    getSolverParams(dfgl).multiproc = false
  end

  # if desired, drawtree in a loop
  treetask, dotreedraw = drawTreeAsyncLoop(tree, opt )

  @info "Do tree based init-inference"
  # if opt.async
  smtasks, hist = taskSolveTreeParametric!(dfgl, tree, oldtree=tree, verbose=verbose, drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs )

  if opt.async && opt.drawtree
    @warn "due to async=true, only keeping task pointer, not stopping the drawtreerate task!  Consider not using .async together with .drawtreerate != 0"
    push!(smtasks, treetask)
  else
    dotreedraw[1] = 0
  end

  @info "Finished tree based Parametric inference"

  return tree, smtasks, hist
end
