
export setVariablePosteriorEstimates!
export updateCliqSolvableDims!, fetchCliqSolvableDims


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

Return true if this clique's down init should be delayed on account of prioritization among sibling separators.

Notes
- process described in issue #344

Dev Notes
- not priorizing order yet (TODO), just avoiding unsolvables at this time.
- Very closely related to getCliqSiblingsPartialNeeds -- refactor likely (NOTE).
- should precompute `allinters`.
"""
function getSiblingsDelayOrder(tree::AbstractBayesTree,
                                cliq::TreeClique,
                                #  prnt,
                                dwnkeys::Vector{Symbol}; # dwinmsgs::LikelihoodMessage;
                                logger=ConsoleLogger())
  # when is a cliq upsolved
  solvedstats = Symbol[:upsolved; :marginalized; :uprecycled]

  # safety net double check
  cliqst = getCliqueStatus(cliq)
  if cliqst in solvedstats
    with_logger(logger) do
      @warn "getSiblingsDelayOrder -- clique status should not be here with a solved cliqst=$cliqst"
    end
    return false
  end

  # get siblings separators
  sibs = getCliqSiblings(tree, cliq, true)
  ids = map(s->s.index, sibs)
  len = length(sibs)
  sibidx = collect(1:len)[ids .== cliq.index][1]
  seps = getCliqSeparatorVarIds.(sibs)
  lielbls = setdiff(ids, cliq.index)
  # get intersect matrix of siblings (should be exactly the same across siblings' csm)
  allinters = Array{Int,2}(undef, len, len)
  dwninters = Vector{Int}(undef, len)
  with_logger(logger) do
    @info "getSiblingsDelayOrder -- number siblings=$(len), sibidx=$sibidx"
  end

  # sum matrix with all "up solved" rows and columns eliminated
  fill!(allinters, 0)
  for i in 1:len
    for j in i:len
      if i != j
        allinters[i,j] = length(intersect(seps[i],seps[j]))
      end
    end
    dwninters[i] = length(intersect(seps[i], dwnkeys))
  end

  # sum "across/over" rows, then columns (i.e. visa versa "along" columns then rows)
  rows = sum(allinters, dims=1)
  cols = sum(allinters, dims=2)

  with_logger(logger) do
      @info "getSiblingsDelayOrder -- allinters=$(allinters)"
      @info "getSiblingsDelayOrder -- rows=$(rows)"
      @info "getSiblingsDelayOrder -- rows=$(cols)"
  end

  # is this clique a non-zero row -- i.e. sum across columns? if not, no further special care needed
  if cols[sibidx] == 0
    with_logger(logger) do
      @info "getSiblingsDelayOrder -- cols[sibidx=$(sibidx))] == 0, no special care needed"
    end
    return false
  end

  # now determine if initializing from below or needdownmsg
  if cliqst in Symbol[:needdownmsg;]
    # be super careful about delay (true) vs pass (false) at this point -- might be partial too TODO
    # return true if delay beneficial to initialization accuracy

    # find which siblings this cliq epends on
    symm = allinters + allinters'
    maskcol = symm[:,sibidx] .> 0
    # lenm = length(maskcol)
    stat = Vector{Symbol}(undef, len)
    stillbusymask = fill(false, len)

    flush(logger.stream)

    # get each sibling status (entering atomic computation segment -- until wait command)
    stat .= getCliqueStatus.(sibs) #[maskcol]

    ## (long down chain case)
    # need different behaviour when all remaining siblings are blocking with :needdownmsg
    remainingmask = stat .== :needdownmsg
    if sum(remainingmask) == length(stat)
      with_logger(logger) do
        @info "getSiblingsDelayOrder -- all blocking: sum(remainingmask) == length(stat), stat=$stat"
      end
      # pick sibling with most overlap in down msgs from parent
      # list of similar length siblings
      candidates = dwninters .== maximum(dwninters)
      if candidates[sibidx]
        # must also pick minimized intersect with other remaing siblings
        maxcan = collect(1:len)[candidates]
        with_logger(logger) do
          @info "getSiblingsDelayOrder -- candidates=$candidates, maxcan=$maxcan, rows=$rows"
        end
        if rows[sibidx] == minimum(rows[maxcan])
          with_logger(logger) do
            @info "getSiblingsDelayOrder -- FORCE DOWN INIT SOLVE ON THIS CLIQUE: $(cliq.index), $(getLabel(cliq))"
          end
          return false
        end
      end
      with_logger(logger) do
        @info "getSiblingsDelayOrder -- not a max and should block"
      end
      return true
    end

    # still busy solving on branches, so potential to delay
    for i in 1:len
      stillbusymask[i] = maskcol[i] && !(stat[i] in solvedstats)
    end
    with_logger(logger) do
        @info "getSiblingsDelayOrder -- busy solving:"
        @info "maskcol=$maskcol"
        @info "stillbusy=$stillbusymask"
    end

    # Too blunt -- should already have returned false by this point perhaps
    if sum(stillbusymask) > 0
      # yes something to delay about
      with_logger(logger) do
        @info "getSiblingsDelayOrder -- yes delay,"
        @info "stat=$stat"
        @info "symm=$symm"
      end
      return true
    end
  end

  with_logger(logger) do
    @info "getSiblingsDelayOrder -- default will not delay"
  end
  flush(logger.stream)
  # carry over default from partial init process
  return false
end


"""
    $SIGNATURES

Return true if both, i.) this clique requires more downward information, ii.) more
downward message information could potentially become available.

Notes
- Delay initialization to the last possible moment.

Dev Notes:

Determine clique truely isn't able to proceed any further:
- should be as self reliant as possible (using clique's status as indicator)
- change status to :mustinitdown if have only partial beliefs so far:
  - combination of status, while partials belief siblings are not :mustinitdown
"""
function getCliqSiblingsPartialNeeds(tree::AbstractBayesTree,
                                      cliq::TreeClique,
                                      #  prnt,
                                      dwinmsgs::LikelihoodMessage;
                                      logger=ConsoleLogger())
  #
  # which incoming messages are partials
  hasPartials = Dict{Symbol, Int}()
  for (sym, tmsg) in dwinmsgs.belief
    # assuming any down message per label that is not partial voids further partial consideration
    if sum(tmsg.inferdim) > 0
      if !haskey(hasPartials, sym)
        hasPartials[sym] = 0
      end
      hasPartials[sym] += 1
    end
  end
  partialKeys = collect(keys(hasPartials))
  
  ## determine who might be able to help init this cliq
  # check sibling separator sets against this clique's separator
  sibs = getCliqSiblings(tree, cliq)
  
  with_logger(logger) do
    @info "getCliqSiblingsPartialNeeds -- CHECK PARTIAL"
  end
  # identify which cliques might have useful information
  localsep = getCliqSeparatorVarIds(cliq)
  seps = Dict{Int, Vector{Symbol}}()
  for si in sibs
    # @show getLabel(si)
    mighthave = intersect(getCliqSeparatorVarIds(si), localsep)
    if length(mighthave) > 0
      seps[si.index] = mighthave
      if getCliqueStatus(si) in [:initialized; :null; :needdownmsg]
        # partials treated special -- this is slightly hacky
        if length(intersect(localsep, partialKeys)) > 0 && length(mighthave) > 0
          # this sibling might have info to delay about
          setCliqDrawColor(cliq,"magenta")
          return true
        end
      end
    end
  end
  # determine if those cliques will / or will not be able to provide more info
  # when does clique change to :mustinitdown
  # default
  return false
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
                                childmsgs=getMsgsUpChildren(csmc, TreeBelief);
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

Related Functions from Upward Inference

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
