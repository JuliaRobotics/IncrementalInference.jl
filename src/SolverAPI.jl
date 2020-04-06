## Various solver API's used in the past.  These functions are due to be standardized, and obsolete code / functions removed.


"""
    $SIGNATURES

Perform inference over the Bayes tree according to `opt::SolverParams`.

Notes
- Variety of options, including fixed-lag solving -- see `getSolverParams(fg)` for details.

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
                    delaycliqs::Vector{Symbol}=Symbol[],
                    recordcliqs::Vector{Symbol}=Symbol[],
                    skipcliqids::Vector{Symbol}=Symbol[],
                    maxparallel::Int=1000,
                    variableOrder::Union{Nothing, Vector{Symbol}}=nothing,
                    variableConstraints::Vector{Symbol}=Symbol[]  ) where G <: DFG.AbstractDFG
  #
  # workaround in case isolated variables occur
  ensureSolvable!(dfgl)
  opt = getSolverParams(dfgl)

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

  orderMethod = 0 < length(variableConstraints) ? :ccolamd : :qr

  # current incremental solver builds a new tree and matches against old tree for recycling.
  tree = wipeBuildNewTree!(dfgl, variableOrder=variableOrder, drawpdf=opt.drawtree, show=opt.showtree, maxparallel=maxparallel,ensureSolvable=false,filepath=joinpath(opt.logpath,"bt.pdf"), variableConstraints=variableConstraints, ordering=orderMethod)
  # setAllSolveFlags!(tree, false)

  # if desired, drawtree in a loop
  treetask, dotreedraw = drawTreeAsyncLoop(tree, opt )

  @info "Do tree based init-inference on tree"
  if opt.async
    smtasks = asyncTreeInferUp!(dfgl, tree, oldtree=oldtree, N=opt.N, drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs )
  else
    smtasks, hist = initInferTreeUp!(dfgl, tree, oldtree=oldtree, N=opt.N, drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs )
  end
  @info "Finished tree based init-inference"

  # NOTE copy of data from new tree in to replace outisde oldtree
  oldtree.bt = tree.bt
  oldtree.btid = tree.btid
  oldtree.cliques = tree.cliques
  oldtree.frontals = tree.frontals
  oldtree.variableOrder = tree.variableOrder
  oldtree.buildTime = tree.buildTime

  if opt.drawtree && opt.async
    @warn "due to async=true, only keeping task pointer, not stopping the drawtreerate task!  Consider not using .async together with .drawtreerate != 0"
    push!(smtasks, treetask)
  else
    dotreedraw[1] = 0
  end

  return oldtree, smtasks, hist
end


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
function solveCliq!(dfgl::G,
                    tree::AbstractBayesTree,
                    cliqid::Symbol;
                    recordcliq::Bool=false,
                    # cliqHistories = Dict{Int,Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}(),
                    maxparallel::Int=50,
                    async::Bool=false  ) where G <: DFG.AbstractDFG
  #
  # hist = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}()
  opt = DFG.getSolverParams(dfgl)

  if opt.isfixedlag
      @info "Quasi fixed-lag is enabled (a feature currently in testing)!"
      fifoFreeze!(dfgl)
  end

  # if !isTreeSolved(treel, skipinitialized=true)
  cliq = whichCliq(tree, cliqid)
  cliqtask = if async
    @async tryCliqStateMachineSolve!(dfgl, tree, cliq.index, drawtree=opt.drawtree, limititers=opt.limititers, downsolve=opt.downsolve,recordcliqs=(recordcliq ? [cliqid] : Symbol[]), incremental=opt.incremental)
  else
    tryCliqStateMachineSolve!(dfgl, tree, cliq.index, drawtree=opt.drawtree, limititers=opt.limititers, downsolve=opt.downsolve,recordcliqs=(recordcliq ? [cliqid] : Symbol[]), incremental=opt.incremental) # N=N
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
                    delaycliqs::Vector{Symbol}=Symbol[],
                    recordcliqs::Vector{Symbol}=Symbol[],
                    skipcliqids::Vector{Symbol}=Symbol[],
                    maxparallel::Int=50)
  #
  @error "Under development, do not use, see #539"
  @info "Solving over the Bayes (Junction) tree."
  smtasks=Vector{Task}()
  hist = Dict{Int, Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}()
  opt = DFG.getSolverParams(dfgl)

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
  smtasks, hist = taskSolveTreeParametric!(dfgl, tree, oldtree=tree, drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs )

  if opt.async && opt.drawtree
    @warn "due to async=true, only keeping task pointer, not stopping the drawtreerate task!  Consider not using .async together with .drawtreerate != 0"
    push!(smtasks, treetask)
  else
    dotreedraw[1] = 0
  end

  @info "Finished tree based Parametric inference"

  return tree, smtasks, hist
end
