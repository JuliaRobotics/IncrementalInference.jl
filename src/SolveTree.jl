
# starting to add exports here


global WORKERPOOL = WorkerPool()

"""
    $SIGNATURES

For use with `multiproc`, nominal use is a worker pool of all processes available above and including 2..., but will return single process [1;] if only the first processes is available.
"""
function setWorkerPool!(pool::Vector{Int}=1 < nprocs() ? setdiff(procs(), [1;]) : [1;])
  global WORKERPOOL
  WORKERPOOL = WorkerPool(pool)
end

function getWorkerPool()
  global WORKERPOOL
  return WORKERPOOL
end

"""
    $SIGNATURES

Previous generation function converting LikelihoodMessage into densities before inference product.

DevNotes
- #913 and consolidation with CSM and `csmc.cliqSubFg`
- FIXME FIXME make sure this is using the proper densities/potentials/likelihoods/factors stored in `csmc.cliqSubFg`.
"""
function packFromIncomingDensities!(dens::Vector{BallTreeDensity},
                                    wfac::Vector{Symbol},
                                    vsym::Symbol,
                                    inmsgs::Array{LikelihoodMessage,1},
                                    manis::T )::Float64 where {T <: Tuple}
  #
  inferdim = 0.0
  for m in inmsgs
    for psym in keys(m.belief)
      if psym == vsym
        pdi = m.belief[vsym]
        push!(dens, manikde!(pdi.val, pdi.bw[:,1], getManifolds(pdi)) )
        push!(wfac, :msg)
        inferdim += pdi.inferdim
      end
      # TODO -- we can inprove speed of search for inner loop
    end
  end

  # return true if at least one of the densities was full dimensional (used for tree based initialization logic)
  return inferdim
end

"""
    $(SIGNATURES)

Add all potentials associated with this clique and vertid to dens.
"""
function packFromLocalPotentials!(dfg::AbstractDFG,
                                  dens::Vector{BallTreeDensity},
                                  wfac::Vector{Symbol},
                                  cliq::TreeClique,
                                  vsym::Symbol,
                                  N::Int,
                                  dbg::Bool=false )
  #
  inferdim = 0.0
  for idfct in getCliqueData(cliq).potentials
    !(exists(dfg, idfct)) && (@warn "$idfct not in clique $(cliq.index)" continue)
    fct = DFG.getFactor(dfg, idfct)
    data = getSolverData(fct)
    # skip partials here, will be caught in packFromLocalPartials!
    if length( findall(getVariableOrder(fct) .== vsym) ) >= 1 && !data.fnc.partial
    # if length( findall(data.fncargvID .== vsym) ) >= 1 && !data.fnc.partial
      p, isinferdim = findRelatedFromPotential(dfg, fct, vsym, N, dbg )
      push!(dens, p)
      push!(wfac, fct.label)
      inferdim += isinferdim
    end
  end

  # return true if at least one of the densities was full dimensional (used for tree based initialization logic)
  return inferdim::Float64
end


function packFromLocalPartials!(fgl::AbstractDFG,
                                partials::Dict{Int, Vector{BallTreeDensity}},
                                cliq::TreeClique,
                                vsym::Symbol,
                                N::Int,
                                dbg::Bool=false  )
  #

  for idfct in getCliqueData(cliq).potentials
    !(exists(fgl, idfct)) && (@warn "$idfct not in clique $(cliq.index)" continue)
    vert = DFG.getFactor(fgl, idfct)
    data = getSolverData(vert)
    if length( findall(getVariableOrder(vert) .== vsym) ) >= 1 && data.fnc.partial
      p, = findRelatedFromPotential(fgl, vert, vsym, N, dbg)
      pardims = data.fnc.usrfnc!.partial
      for dimnum in pardims
        if haskey(partials, dimnum)
          push!(partials[dimnum], marginal(p,[dimnum]))
        else
          partials[dimnum] = BallTreeDensity[marginal(p,[dimnum])]
        end
      end
    end
  end
  nothing
end



"""
    $(SIGNATURES)

Perform one step of the minibatch clique Gibbs operation for solving the Chapman-Kolmogov
trasit integral -- here involving separate approximate functional convolution and
product operations.
"""
function cliqGibbs( dfg::AbstractDFG,
                    cliq::TreeClique,
                    vsym::Symbol,
                    inmsgs::Array{LikelihoodMessage,1},
                    N::Int,
                    dbg::Bool,
                    manis::Tuple,
                    logger=ConsoleLogger()  )
  #
  # several optimizations can be performed in this function TODO

  # consolidate NBPMessages and potentials
  dens = Array{BallTreeDensity,1}()
  partials = Dict{Int, Vector{BallTreeDensity}}()
  wfac = Vector{Symbol}()

  inferdim = 0.0
  inferdim += packFromIncomingDensities!(dens, wfac, vsym, inmsgs, manis)
  inferdim += packFromLocalPotentials!(dfg, dens, wfac, cliq, vsym, N)
  packFromLocalPartials!(dfg, partials, cliq, vsym, N, dbg)

  potprod = !dbg ? nothing : PotProd(vsym, getVal(dfg,vsym), Array{Float64,2}(undef, 0,0), dens, wfac)
      # pts,inferdim = predictbelief(dfg, vsym, useinitfct)  # for reference only
  pGM = productbelief(dfg, vsym, dens, partials, N, dbg=dbg, logger=logger )
  if dbg  potprod.product = pGM  end


  # @info " "
  return pGM, potprod, inferdim
end

"""
    $SIGNATURES

Dev Notes
- part of refactoring fmcmc.
- function seems excessive
"""
function compileFMCMessages(fgl::AbstractDFG,
                            lbls::Vector{Symbol},
                            logger=ConsoleLogger())
  #
  d = Dict{Symbol,TreeBelief}()
  for vsym in lbls
    vari = DFG.getVariable(fgl,vsym)
    pden = getKDE(vari)
    bws = vec(getBW(pden)[:,1])
    manis = getSofttype(vari) |> getManifolds
    d[vsym] = TreeBelief(vari) # getVal(vari), bws, manis, getSolverData(vari).inferdim
    with_logger(logger) do
      @info "fmcmc! -- getSolverData(vari=$(vari.label)).inferdim=$(getSolverData(vari).inferdim)"
    end
  end
  return d
end

# global countsolve = 0

function doFMCIteration(fgl::AbstractDFG,
                        vsym::Symbol,
                        cliq::TreeClique,
                        fmsgs,
                        N::Int,
                        dbg::Bool,
                        logger=ConsoleLogger()  )
  #
# global countsolve
  vert = DFG.getVariable(fgl, vsym)
  if !getSolverData(vert).ismargin
    # we'd like to do this more pre-emptive and then just execute -- just point and skip up only msgs
    densPts, potprod, inferdim = cliqGibbs(fgl, cliq, vsym, fmsgs, N, dbg, getSofttype(vert) |> getManifolds, logger)

      # countsolve += 1
      # saveDFG(fgl, "/tmp/fix/$(vsym)_$(countsolve)")

    if size(densPts,1)>0
      updvert = DFG.getVariable(fgl, vsym)  # TODO --  can we remove this duplicate getVert?
      setValKDE!(updvert, densPts, true, inferdim)
      # FIXME Restore debugging inside mm solve
      # if dbg
      #   push!(dbgvals.prods, potprod)
      #   push!(dbgvals.lbls, Symbol(updvert.label))
      # end
    end
  end
  nothing
end

"""
    $(SIGNATURES)

Iterate successive approximations of clique marginal beliefs by means
of the stipulated proposal convolutions and products of the functional objects
for tree clique `cliq`.
"""
function fmcmc!(fgl::AbstractDFG,
                cliq::TreeClique,
                fmsgs::Vector{LikelihoodMessage},
                lbls::Vector{Symbol},
                N::Int,
                MCMCIter::Int,
                dbg::Bool=false,
                logger=ConsoleLogger(),
                multithreaded::Bool=false  )
  #
  with_logger(logger) do
    @info "---------- successive fnc approx ------------$(getLabel(cliq))"
  end
  # repeat several iterations of functional Gibbs sampling for fixed point convergence
  if length(lbls) == 1
      MCMCIter=1
  end
  mcmcdbg = Array{CliqGibbsMC,1}()

  for iter in 1:MCMCIter
    # iterate through each of the variables, KL-divergence tolerence would be nice test here
    with_logger(logger) do
      @info "#$(iter)\t -- "
    end
    dbgvals = !dbg ? nothing : CliqGibbsMC([], Symbol[])

    for vsym in lbls
        doFMCIteration(fgl, vsym, cliq, fmsgs, N, dbg, logger)
    end
    !dbg ? nothing : push!(mcmcdbg, dbgvals)
  end

  # populate dictionary for return NBPMessage in multiple dispatch
  msgdict = compileFMCMessages(fgl, lbls, logger)

  return mcmcdbg, msgdict
end


"""
    $(SIGNATURES)

Perform computations required for the upward message passing during belief propation on the Bayes (Junction) tree.
This function is usually called as via remote_call for multiprocess dispatch.

Example
```julia
inp = ExploreTreeType(fg,tree,cliq,parent,childmsgs)
urt = upGibbsCliqueDensity(inp)
```
- `fg` factor graph,
- `tree` Bayes tree,
- `cliq` which cliq to perform the computation on,
- `parent` the parent clique to where the upward message will be sent,
- `childmsgs` is for any incoming messages from child cliques.
"""
function upGibbsCliqueDensity(inp::FullExploreTreeType{T,T2},
                              N::Int=100,
                              dbg::Bool=false,
                              iters::Int=3,
                              logger=ConsoleLogger()  ) where {T, T2}
  #
  with_logger(logger) do
    @info "up w $(length(inp.sendmsgs)) msgs"
  end
  # Local mcmc over belief functions
  # this is so slow! TODO Can be ignored once we have partial working
  # loclfg = nprocs() < 2 ? deepcopy(inp.fg) : inp.fg

  # TODO -- some weirdness with: d,. = d = ., nothing
  mcmcdbg = Array{CliqGibbsMC,1}()
  d = Dict{Symbol,TreeBelief}()

  priorprods = Vector{CliqGibbsMC}()

  cliqdata = getCliqueData(inp.cliq)

  with_logger(logger) do
    for el in inp.sendmsgs, (id,msg) in el.belief
      @info "inp.sendmsgs[$id].inferdim=$(msg.inferdim)"
    end
  end

  # use nested structure for more fficient Chapman-Kolmogorov solution approximation
  if false
    IDS = [cliqdata.frontalIDs;cliqdata.separatorIDs] #inp.cliq.attributes["frontalIDs"]
    mcmcdbg, d = fmcmc!(inp.fg, inp.cliq, inp.sendmsgs, IDS, N, iters, dbg, logger)
  else
    # NOTE -- previous mistake, must iterate over directsvarIDs also (or incorporate once at the right time)
    dummy, d = fmcmc!(inp.fg, inp.cliq, inp.sendmsgs, cliqdata.directFrtlMsgIDs, N, 1, dbg, logger, true)
    if length(cliqdata.msgskipIDs) > 0
      dummy, dd = fmcmc!(inp.fg, inp.cliq, inp.sendmsgs, cliqdata.msgskipIDs, N, 1, dbg, logger, true)
      for md in dd d[md[1]] = md[2]; end
    end
    if length(cliqdata.itervarIDs) > 0
      mcmcdbg, ddd = fmcmc!(inp.fg, inp.cliq, inp.sendmsgs, cliqdata.itervarIDs, N, iters, dbg, logger, false)
      for md in ddd d[md[1]] = md[2]; end
    end
    if length(cliqdata.directPriorMsgIDs) > 0
      doids = setdiff(cliqdata.directPriorMsgIDs, cliqdata.msgskipIDs)
      priorprods, dddd = fmcmc!(inp.fg, inp.cliq, inp.sendmsgs, doids, N, 1, dbg, logger, true)
      for md in dddd d[md[1]] = md[2]; end
    end
  end

  return d
end



"""
    $SIGNATURES

Perform Chapman-Kolmogorov transit integral approximation for `cliq` in downward pass direction.

Notes
- Only update frontal variables of the clique.
"""
function downGibbsCliqueDensity(fg::AbstractDFG,
                                cliq::TreeClique,
                                dwnMsgs::Array{LikelihoodMessage,1},
                                N::Int=100,
                                MCMCIter::Int=3,
                                dbg::Bool=false,
                                usemsgpriors::Bool=false,
                                logger=ConsoleLogger() )
  #
  # TODO standardize function call to have similar stride to upGibbsCliqueDensity
  # @info "down"
  with_logger(logger) do
    @info "cliq=$(cliq.index), downGibbsCliqueDensity -- going for down fmcmc"
  end
  fmmsgs = usemsgpriors ? Array{LikelihoodMessage,1}() : dwnMsgs
  frtls = getFrontals(cliq)

  # TODO, do better check if there is structure between multiple frontals
  niters = length(frtls) == 1 ? 1 : MCMCIter
  # TODO standize with upsolve and variable solver order
  mcmcdbg, d = fmcmc!(fg, cliq, fmmsgs, frtls, N, niters, dbg)
  m = dwnPrepOutMsg(fg, cliq, dwnMsgs, d, logger)

  outmsglbl = Dict{Symbol, Int}()
  if dbg
    for (ke, va) in m.belief
      outmsglbl[Symbol(fg.g.vertices[ke].label)] = ke
    end
  end

  with_logger(logger) do
    @info "cliq=$(cliq.index), downGibbsCliqueDensity -- convert to BallTreeDensities."
  end

  # Always keep dwn messages in cliq data
  dwnkeepmsgs = LikelihoodMessage()
  for (msgsym, val) in m.belief
    dwnkeepmsgs.belief[msgsym] = convert(Tuple{BallTreeDensity,Float64}, val)
  end
  setDwnMsg!(cliq, dwnkeepmsgs)

  # down solving complete, set flag
  getCliqueData(cliq).downsolved = true

  mdbg = !dbg ? DebugCliqMCMC() : DebugCliqMCMC(mcmcdbg, m, outmsglbl, CliqGibbsMC[])
  with_logger(logger) do
    @info "cliq=$(cliq.index), downGibbsCliqueDensity -- finished."
  end
  return DownReturnBPType(m, mdbg, d, dwnkeepmsgs)
end
function downGibbsCliqueDensity(fg::AbstractDFG,
                                cliq::TreeClique,
                                dwnMsgs::LikelihoodMessage,
                                N::Int=100,
                                MCMCIter::Int=3,
                                dbg::Bool=false,
                                usemsgpriors::Bool=false,
                                logger=ConsoleLogger() )
  #
  with_logger(logger) do
    @info "cliq=$(cliq.index), downGibbsCliqueDensity -- convert BallTreeDensities to LikelihoodMessage."
  end
  ind = Dict{Symbol, TreeBelief}()
  sflbls = listVariables(fg)
  for (lbl, bel) in dwnMsgs.belief
    if lbl in sflbls
      ind[lbl] = TreeBelief(bel[1], bel[2], getSofttype(getVariable(fg, lbl)))
    end
  end
  ndms = LikelihoodMessage[LikelihoodMessage(ind);]
  with_logger(logger) do
    @info "cliq=$(cliq.index), downGibbsCliqueDensity -- call with LikelihoodMessage."
  end
  downGibbsCliqueDensity(fg, cliq, ndms, N, MCMCIter, dbg, usemsgpriors, logger)
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
  fg_ = csmc.cliqSubFg
  tree_ = csmc.tree
  cliq = csmc.cliq

  # ?? TODO use subgraph copy of factor graph for operations and transfer frontal variables only

  with_logger(logger) do
    @info "=== start Clique $(getLabel(cliq)) ======================"
  end
  ett = FullExploreTreeType(fg_, nothing, cliq, nothing, childmsgs)

  if multiproc
    cliqc = deepcopy(cliq)
    btnd = getCliqueData(cliqc)
    # redirect to new unused so that CAN be serialized
    btnd.upMsgChannel = Channel{LikelihoodMessage}(1)
    # btnd.initDownChannel = Channel{LikelihoodMessage}(1) # TODO deprecate
    btnd.dwnMsgChannel = Channel{LikelihoodMessage}(1)
    btnd.solveCondition = Condition()
    ett.cliq = cliqc
    # TODO create new dedicate file for separate process to log with
    try
      retdict = remotecall_fetch(upGibbsCliqueDensity, getWorkerPool(), ett, N, dbg, iters)
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
    retdict = upGibbsCliqueDensity(ett, N, dbg, iters, logger)
  end

  with_logger(logger) do
    @info "=== end Clique $(getLabel(cliq)) ========================"
  end

  return retdict
end



"""
    $SIGNATURES

Special internal function to try return the clique data if succesfully identified in `othertree::AbstractBayesTree`,
based on contents of `seeksSimilar::BayesTreeNodeData`.

Notes
- Used to identify and skip similar cliques (i.e. recycle computations)
"""
function attemptTreeSimilarClique(othertree::AbstractBayesTree, 
                                  seeksSimilar::BayesTreeNodeData  )
  #
  # inner convenience function for returning empty clique
  function EMPTYCLIQ()
    clq = TreeClique(-1,"null")
    setLabel!(clq, "")
    setCliqueData!(clq, BayesTreeNodeData())
    return clq
  end

  # does the other clique even exist?
  seekFrontals = getCliqFrontalVarIds(seeksSimilar)
  if !hasClique(othertree, seekFrontals[1])
    return EMPTYCLIQ()
  end

  # do the cliques share the same frontals?
  otherCliq = getClique(othertree, seekFrontals[1])
  otherFrontals = getCliqFrontalVarIds(otherCliq)
  commonFrontals = intersect(seekFrontals, otherFrontals)
  if length(commonFrontals) != length(seekFrontals) || length(commonFrontals) != length(otherFrontals)
    return EMPTYCLIQ()
  end

  # do the cliques share the same separator variables?
  seekSeparator = getCliqSeparatorVarIds(seeksSimilar)
  otherSeparator = getCliqSeparatorVarIds(otherCliq)
  commonSep = intersect(seekSeparator, otherSeparator)
  if length(commonSep) != length(seekSeparator) || length(commonSep) != length(otherSeparator)
    return EMPTYCLIQ()
  end

  # do the cliques use the same factors (potentials)
  seekPotentials = getCliqFactorIds(seeksSimilar)
  otherFactors = getCliqFactorIds(otherCliq)
  commonFactors = intersect(seekPotentials, otherFactors)
  if length(commonFactors) != length(seekPotentials) || length(commonFactors) != length(otherFactors)
    return EMPTYCLIQ()
  end

  # lets assume they are the same
  return otherCliq::TreeClique
end




#