
export setVariablePosteriorEstimates!
export updateCliqSolvableDims!, fetchCliqSolvableDims

# using Serialization

"""
    $SIGNATURES

Internal helper function to save a dfg object to LogPath during clique state machine operations.

Notes
- will only save dfg object if `opts.dbg=true`

Related

saveDFG, loadDFG!, loadDFG
"""
function _dbgCSMSaveSubFG(csmc::CliqStateMachineContainer, filename::String)

  opt = getSolverParams(csmc.cliqSubFg)
  
  if opt.dbg
    folder::String=joinpath(opt.logpath,"logs","cliq$(csmc.cliq.index)")
    if !ispath(folder)
      mkpath(folder)
    end
    # NOTE there was a bug using saveDFG, so used serialize, left for future use  
    # serialize(joinpath(folder, filename), csmc.cliqSubFg)
    DFG.saveDFG(csmc.cliqSubFg, joinpath(folder, filename))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(folder, "$(filename).pdf"))
  end
  
  opt.dbg

end


"""
    $SIGNATURES

Set all up `upsolved` and `downsolved` cliq data flags `to::Bool=false`.
"""
function setAllSolveFlags!(treel::AbstractBayesTree, to::Bool=false)::Nothing
  for (id, cliq) in getCliques(treel)
    cliqdata = getCliqueData(cliq)
    setCliqueStatus!(cliqdata, :null)
    cliqdata.upsolved = to
    cliqdata.downsolved = to
  end
  nothing
end

"""
    $SIGNATURES

Return true or false depending on whether the tree has been fully initialized/solved/marginalized.
"""
function isTreeSolved(treel::AbstractBayesTree; skipinitialized::Bool=false)
  acclist = Symbol[:upsolved; :downsolved; :marginalized]
  skipinitialized ? nothing : push!(acclist, :initialized)
  for (clid, cliq) in getCliques(treel)
    if !(getCliqueStatus(cliq) in acclist)
      return false
    end
  end
  return true
end

function isTreeSolvedUp(treel::AbstractBayesTree)
  for (clid, cliq) in getCliques(treel)
    if getCliqueStatus(cliq) != :upsolved
      return false
    end
  end
  return true
end


"""
    $SIGNATURES

Return `::Bool` on whether all variables in this `cliq` are marginalzed.
"""
function isCliqMarginalizedFromVars(subfg::AbstractDFG, cliq::TreeClique)
  for vert in getCliqVars(subfg, cliq)
    if !isMarginalized(vert)
      return false
    end
  end
  return true
end


"""
    $SIGNATURES

Reset the Bayes (Junction) tree so that a new upsolve can be performed.

Notes
- Will change previous clique status from `:downsolved` to `:initialized` only.
- Sets the color of tree clique to `lightgreen`.
"""
function resetTreeCliquesForUpSolve!(treel::AbstractBayesTree)::Nothing
  acclist = Symbol[:downsolved;]
  for (clid, cliq) in getCliques(treel)
    if getCliqueStatus(cliq) in acclist
      setCliqueStatus!(cliq, :initialized)
      setCliqDrawColor(cliq, "sienna")
    end
  end
  nothing
end


"""
    $SIGNATURES

Store/cache a clique's solvable dimensions.
"""
function updateCliqSolvableDims!( cliq::TreeClique,
                                  sdims::Dict{Symbol, Float64},
                                  logger=ConsoleLogger() )::Nothing
  #
  cliqd = getCliqueData(cliq)
  csd = getSolvableDims(cliqd)
  if isready(csd)
    take!(csd)
    with_logger(logger) do
      @info "cliq $(cliq.index), updateCliqSolvableDims! -- cleared solvableDims"
    end
  end
  put!(csd, sdims)
  with_logger(logger) do
      @info "cliq $(cliq.index), updateCliqSolvableDims! -- updated"
  end
  nothing
end


"""
    $SIGNATURES

Retrieve a clique's cached solvable dimensions (since last update).
"""
function fetchCliqSolvableDims(cliq::TreeClique)::Dict{Symbol,Float64}
  cliqd = getCliqueData(cliq)
  csd = getSolvableDims(cliqd)
  if isready(csd)
    return csd.data[1]
  end
  return fetch(csd)
end


"""
    $SIGNATURES

Return true there is no other sibling that will make progress.

Notes
- Relies on sibling priority order with only one "currently best" option that will force progress in global upward inference.
- Return false if one of the siblings is still busy
"""
function areSiblingsRemaingNeedDownOnly(tree::AbstractBayesTree,
                                        cliq::TreeClique  )::Bool
  #
  stillbusylist = [:null; :initialized;]
  prnt = getParent(tree, cliq)
  if length(prnt) > 0
    for si in getChildren(tree, prnt[1])
      # are any of the other siblings still busy?
      if si.index != cliq.index && getCliqueStatus(si) in stillbusylist
        return false
      end
    end
  end

  # nope, everybody is waiting for something to change -- proceed with forcing a cliq solve
  return true
end


"""
    $SIGNATURES

Approximate Chapman-Kolmogorov transit integral and return separator marginals as messages to pass up the Bayes (Junction) tree, along with additional clique operation values for debugging.

Notes
- `onduplicate=true` by default internally uses deepcopy of factor graph and Bayes tree, and does **not** update the given objects.  Set false to update `fgl` and `treel` during compute.

Future
- TODO: internal function chain is too long and needs to be refactored for maintainability.
"""
function approxCliqMarginalUp!( csmc::CliqStateMachineContainer,
                                childmsgs=fetchMsgsUpChildren(csmc, TreeBelief);
                                N::Int=getSolverParams(csmc.cliqSubFg).N,
                                dbg::Bool=getSolverParams(csmc.cliqSubFg).dbg,
                                multiproc::Bool=getSolverParams(csmc.cliqSubFg).multiproc,
                                logger=ConsoleLogger(),
                                iters::Int=3,
                                drawpdf::Bool=false  )
  #
  # use subgraph copy of factor graph for operations and transfer variables results later only
  fg_ = csmc.cliqSubFg
  tree_ = csmc.tree
  cliq = csmc.cliq

  with_logger(logger) do
    @info "=== start Clique $(getLabel(cliq)) ======================"
  end

  if multiproc
    cliqc = deepcopy(cliq)
    btnd = getCliqueData(cliqc)
    # redirect to new unused so that CAN be serialized
    btnd.upMsgChannel = Channel{LikelihoodMessage}(1)
    btnd.dwnMsgChannel = Channel{LikelihoodMessage}(1)
    btnd.solveCondition = Condition()
    # ett.cliq = cliqc
    # TODO create new dedicated file for separate process to log with
    try
      retdict = remotecall_fetch(upGibbsCliqueDensity, getWorkerPool(), fg_, cliqc, childmsgs, N, dbg, iters)
    catch ex
      with_logger(logger) do
        @info ex
        @error ex
        flush(logger.stream)
        msg = sprint(showerror, ex)
        @error msg
      end
      flush(logger.stream)
      error(ex)
    end
  else
    with_logger(logger) do
      @info "Single process upsolve clique=$(cliq.index)"
    end
    retdict = upGibbsCliqueDensity(fg_, cliq, childmsgs, N, dbg, iters, logger)
  end

  with_logger(logger) do
    @info "=== end Clique $(getLabel(cliq)) ========================"
  end

  return retdict
end


"""
    $SIGNATURES

Determine which variables to iterate or compute directly for downward tree pass of inference.

DevNotes
- # TODO see #925

Related

directPriorMsgIDs, directFrtlMsgIDs, directAssignmentIDs, mcmcIterationIDs
"""
function determineCliqVariableDownSequence( subfg::AbstractDFG, 
                                            cliq::TreeClique; 
                                            solvable::Int=1, 
                                            logger=ConsoleLogger())
  #
  frtl = getCliqFrontalVarIds(cliq)

  adj, varLabels, FactorLabels = DFG.getBiadjacencyMatrix(subfg, solvable=solvable)
  mask = (x-> x in frtl).(varLabels)
  newFrtlOrder = varLabels[mask]
  subAdj = adj[:,mask]
    #TODO don't use this getAdjacencyMatrixSymbols, #604
    # adj = DFG.getAdjacencyMatrixSymbols(subfg, solvable=solvable)
    # mask = map(x->(x in frtl), adj[1,:])
    # subAdj = adj[2:end,mask] .!= nothing
    # newFrtlOrder = Symbol.(adj[1,mask])

  crossCheck = 1 .< sum(Int.(subAdj), dims=2)
  iterVars = Symbol[]
  for i in 1:length(crossCheck)
    # must add associated variables to iterVars
    if crossCheck[i]
          # # DEBUG loggin
          # with_logger(logger) do
          #   @info "newFrtlOrder=$newFrtlOrder"
          #   @info "(subAdj[i,:]).nzind=$((subAdj[i,:]).nzind)"
          # end
          # flush(logger.stream)
      # find which variables are associated
      varSym = newFrtlOrder[(subAdj[i,:]).nzind]
      union!(iterVars, varSym)
    end
  end

  # return iteration list ordered by frtl
  return intersect(frtl, iterVars)
end


"""
    $SIGNATURES

Perform downward direction solves on a sub graph fragment.
Calculates belief on each of the frontal variables and iterate if required.

Notes
- uses all factors connected to the frontal variables.
- assumes `subfg` was properly prepared before calling.
- has multi-process option.

Dev Notes
- TODO incorporate variation possible due to cross frontal factors.
- cleanup and updates required, and @spawn jl 1.3
"""
function solveCliqDownFrontalProducts!( subfg::AbstractDFG,
                                        cliq::TreeClique,
                                        opts::SolverParams,
                                        logger=ConsoleLogger();
                                        MCIters::Int=3 )
  #
  # get frontal variables for this clique
  frsyms = getCliqFrontalVarIds(cliq)

  # determine if cliq has cross frontal factors
  # iterdwn, directdwns, passmsgs?
  iterFrtls = determineCliqVariableDownSequence(subfg,cliq, logger=logger)

  # direct frontals
  directs = setdiff(frsyms, iterFrtls)

  # ignore limited fixed lag variables
  fixd = map(x->opts.limitfixeddown && isMarginalized(subfg,x), frsyms)
  skip = frsyms[fixd]
  iterFrtls = setdiff(iterFrtls, skip)
  directs = setdiff(directs, skip)
  with_logger(logger) do
    @info "cliq $(cliq.index), solveCliqDownFrontalProducts!, skipping marginalized keys=$(skip)"
  end


  # use new localproduct approach
  if opts.multiproc
    downresult = Dict{Symbol, Tuple{BallTreeDensity, Float64, Vector{Symbol}}}()
    @sync for i in 1:length(directs)
      @async begin
        downresult[directs[i]] = remotecall_fetch(localProductAndUpdate!, getWorkerPool(), subfg, directs[i], false)
        # downresult[directs[i]] = remotecall_fetch(localProductAndUpdate!, upp2(), subfg, directs[i], false)
      end
    end
    with_logger(logger) do
      @info "cliq $(cliq.index), solveCliqDownFrontalProducts!, multiproc keys=$(keys(downresult))"
    end
    for fr in directs
        with_logger(logger) do
            @info "cliq $(cliq.index), solveCliqDownFrontalProducts!, key=$(fr), infdim=$(downresult[fr][2]), lbls=$(downresult[fr][3])"
        end
      setValKDE!(subfg, fr, downresult[fr][1], false, downresult[fr][2])
    end
    for mc in 1:MCIters, fr in iterFrtls
      try
        result = remotecall_fetch(localProductAndUpdate!, getWorkerPool(), subfg, fr, false)
        # result = remotecall_fetch(localProductAndUpdate!, upp2(), subfg, fr, false)
        setValKDE!(subfg, fr, result[1], false, result[2])
        with_logger(logger) do
          @info "cliq $(cliq.index), solveCliqDownFrontalProducts!, iter key=$(fr), infdim=$(result[2]), lbls=$(result[3])"
        end
      catch ex
        # what if results contains an error?
        with_logger(logger) do
          @error ex
          flush(logger.stream)
          msg = sprint(showerror, ex)
          @error msg
        end
        error(ex)
      end
    end
  else
    # do directs first
    for fr in directs
      localProductAndUpdate!(subfg, fr, true, logger)
    end
    #do iters next
    for mc in 1:MCIters, fr in iterFrtls
      localProductAndUpdate!(subfg, fr, true, logger)
    end
  end

  return nothing
end





#
