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
                    oldtree::BayesTree=emptyBayesTree();
                    delaycliqs::Vector{Symbol}=Symbol[],
                    recordcliqs::Vector{Symbol}=Symbol[],
                    skipcliqids::Vector{Symbol}=Symbol[],
                    maxparallel::Int=50  ) where G <: DFG.AbstractDFG
  #
  @info "Solving over the Bayes (Junction) tree."
  smtasks=Vector{Task}()
  hist = Dict{Int, Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}()
  opt = DFG.getSolverParams(dfgl)

  if opt.isfixedlag
      @info "Quasi fixed-lag is enabled (a feature currently in testing)!"
      fifoFreeze!(dfgl)
  end

  # current incremental solver builds a new tree and matches against old tree for recycling.
  tree = wipeBuildNewTree!(dfgl, drawpdf=opt.drawtree, show=opt.showtree, maxparallel=maxparallel, filepath=joinpath(getSolverParams(dfgl).logpath,"bt.pdf"))
  # setAllSolveFlags!(tree, false)

  @info "Do tree based init-inference on tree"
  if opt.async
    smtasks = asyncTreeInferUp!(dfgl, tree, oldtree=oldtree, N=opt.N, drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs )
  else
    smtasks, hist = initInferTreeUp!(dfgl, tree, oldtree=oldtree, N=opt.N, drawtree=opt.drawtree, recordcliqs=recordcliqs, limititers=opt.limititers, downsolve=opt.downsolve, incremental=opt.incremental, skipcliqids=skipcliqids, delaycliqs=delaycliqs )
  end
  @info "Finished tree based init-inference"

  # transfer new tree to outside parameter
  oldtree.bt = tree.bt
  oldtree.btid = tree.btid
  oldtree.cliques = tree.cliques
  oldtree.frontals = tree.frontals

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
                    tree::BayesTree,
                    cliqid::Symbol;
                    recordcliq::Bool=false,
                    # cliqHistories = Dict{Int,Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}(),
                    maxparallel::Int=50  ) where G <: DFG.AbstractDFG
  #
  # hist = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}()
  opt = DFG.getSolverParams(dfgl)

  if opt.isfixedlag
      @info "Quasi fixed-lag is enabled (a feature currently in testing)!"
      fifoFreeze!(dfgl)
  end

  # if !isTreeSolved(treel, skipinitialized=true)
  cliq = whichCliq(tree, cliqid)
  cliqtask = @async tryCliqStateMachineSolve!(dfgl, tree, cliq.index, drawtree=opt.drawtree, limititers=opt.limititers, downsolve=opt.downsolve,recordcliqs=(recordcliq ? [cliqid] : Symbol[]), incremental=opt.incremental) # N=N
  # end # if


  # post-hoc store possible state machine history in clique (without recursively saving earlier history inside state history)
  # assignTreeHistory!(tree, cliqHistories)

  # cliqHistories
  return cliqtask
end



"""
    $SIGNATURES

Perform up and down message passing (multi-process) algorithm for full sum-product solution of all continuous marginal beliefs.

Notes
- For legacy versions of tree traversal, see `inferOverTreeIterative!` instead.
"""
function inferOverTree!(dfg::G,
                        bt::BayesTree;
                        oldtree::BayesTree=emptyBayesTree(),
                        N::Int=100,
                        upsolve::Bool=true,
                        downsolve::Bool=true,
                        dbg::Bool=false,
                        drawpdf::Bool=false,
                        treeinit::Bool=false,
                        incremental::Bool=false,
                        limititers::Int=1000,
                        skipcliqids::Vector{Symbol}=Symbol[],
                        recordcliqs::Vector{Symbol}=Symbol[]  ) where G <: AbstractDFG
  #

  @info "Solving over the Bayes (Junction) tree."
  smtasks=Vector{Task}()
  ch = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}()
  setAllSolveFlags!(bt, false)

  @info "Do tree based init-inference on tree"
  if dbg
    smtasks = asyncTreeInferUp!(dfg, bt, oldtree=oldtree, N=N, drawtree=drawpdf, recordcliqs=recordcliqs, limititers=limititers, downsolve=downsolve, incremental=incremental, skipcliqids=skipcliqids )
  else
    smtasks, ch = initInferTreeUp!(dfg, bt, oldtree=oldtree, N=N, drawtree=drawpdf, recordcliqs=recordcliqs, limititers=limititers, downsolve=downsolve, incremental=incremental, skipcliqids=skipcliqids )
  end
  @info "Finished tree based init-inference"

  return smtasks, ch
end

"""
    $SIGNATURES

Perform up and down message passing (multi-process) algorithm for full sum-product solution of all continuous marginal beliefs.

Notes
- Legacy support function, use `inferOverTree!` instead that is based on the state machine method.
- Previous versions of the code used iterative loops to traverse the Bayes (Junction) tree.
- Even older code is available as `inferOverTreeR!`
"""
function inferOverTreeIterative!(dfg::G,
                                 bt::BayesTree;
                                 N::Int=100,
                                 dbg::Bool=false,
                                 drawpdf::Bool=false  ) where G <: AbstractDFG
  #
  # @info "Batch rather than incremental solving over the Bayes (Junction) tree."
  # setAllSolveFlags!(bt, false)
  @info "Ensure all nodes are initialized"
  ensureAllInitialized!(dfg)
  @info "Do multi-process upward pass of inference on tree"
  upMsgPassingIterative!(ExploreTreeType(dfg, bt, bt.cliques[1], nothing, NBPMessage[]),N=N, dbg=dbg, drawpdf=drawpdf);
  @info "Do multi-process downward pass of inference on tree"
  downMsgPassingIterative!(ExploreTreeType(dfg, bt, bt.cliques[1], nothing, NBPMessage[]),N=N, dbg=dbg, drawpdf=drawpdf);
  return smtasks, ch
end

"""
    $SIGNATURES

Perform up and down message passing (single process, recursive) algorithm for full sum-product solution of all continuous marginal beliefs.
"""
function inferOverTreeR!(fgl::G,
                         bt::BayesTree;
                         N::Int=100,
                         dbg::Bool=false,
                         drawpdf::Bool=false,
                         treeinit::Bool=false  )::Tuple{Vector{Task},Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}} where G <: AbstractDFG
  #
  @info "Batch rather than incremental solving over the Bayes (Junction) tree."
  setAllSolveFlags!(bt, false)
  @info "Ensure all nodes are initialized"
  smtasks = Vector{Task}()
  ch = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}()
  if treeinit
    smtasks, ch = initInferTreeUp!(fgl, bt, N=N, drawtree=drawpdf)
  else
    @info "Do conventional recursive up inference over tree"
    ensureAllInitialized!(fgl)
    upMsgPassingRecursive(ExploreTreeType(fgl, bt, bt.cliques[1], nothing, NBPMessage[]), N=N, dbg=dbg, drawpdf=drawpdf);
  end
  @info "Do recursive down inference over tree"
  downMsgPassingRecursive(ExploreTreeType(fgl, bt, bt.cliques[1], nothing, NBPMessage[]), N=N, dbg=dbg, drawpdf=drawpdf);
  return smtasks, ch
end





"""
    $(SIGNATURES)

Perform multimodal incremental smoothing and mapping (mm-iSAM) computations over given factor graph `fgl::FactorGraph` on the local computer.  A pdf of the Bayes (Junction) tree will be generated in the working folder with `drawpdf=true`
"""
function batchSolve!(dfg::G,
                     oldtree::BayesTree=emptyBayesTree();
                     upsolve::Bool=true,
                     downsolve::Bool=true,
                     drawpdf::Bool=false,
                     show::Bool=false,
                     N::Int=100,
                     recursive::Bool=false,
                     dbg::Bool=false,
                     treeinit::Bool=false,
                     incremental::Bool=false,
                     limititers::Int=1000,
                     skipcliqids::Vector{Symbol}=Symbol[],
                     recordcliqs::Vector{Symbol}=Symbol[],
                     returntasks::Bool=false,
                     maxparallel::Int=50  ) where G <: AbstractDFG
  #
  @warn "deprecated batchSolve! in favor of new solveTree! interface with the same and more functionality."
  if DFG.getSolverParams(dfg).isfixedlag
      @info "Quasi fixed-lag is enabled (a feature currently in testing)!"
      fifoFreeze!(dfg)
  end
  tree = wipeBuildNewTree!(dfg, drawpdf=false, maxparallel=maxparallel)
  drawpdf ? drawTree(tree, show=show) : nothing
  # show ? showTree() : nothing

  smtasks = Vector{Task}()
  if recursive
    # recursive is a single core method that is slower but occasionally helpful for better stack traces during debugging
    smtasks, ch = inferOverTreeR!(dfg, tree, N=N, drawpdf=drawpdf, dbg=dbg, treeinit=treeinit)
  else
    smtasks, ch = inferOverTree!(dfg, tree, oldtree=oldtree, N=N, drawpdf=drawpdf, dbg=dbg, treeinit=treeinit,
                                 limititers=limititers, recordcliqs=recordcliqs, upsolve=upsolve,
                                 downsolve=downsolve, incremental=incremental, skipcliqids=skipcliqids  )
  end

  # later development allows tasks for each cliq state machine to be returned also
  return tree, smtasks
end
