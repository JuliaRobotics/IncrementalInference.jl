
# starting to add exports here
export findRelatedFromPotential


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
    $(SIGNATURES)

Compute proposal belief on `vertid` through `fct` representing some constraint in factor graph.
Always full dimension variable node -- partial constraints will only influence subset of variable dimensions.
The remaining dimensions will keep pre-existing variable values.

Notes
- fulldim is true when "rank-deficient" -- TODO swap to false (or even float)
"""
function findRelatedFromPotential(dfg::AbstractDFG,
                                  fct::DFGFactor,
                                  varid::Symbol,
                                  N::Int,
                                  dbg::Bool=false;
                                  solveKey::Symbol=:default )
  #
  # assuming it is properly initialized TODO
  ptsbw = evalFactor2(dfg, fct, varid, solveKey=solveKey, N=N, dbg=dbg);
  # determine if evaluation is "dimension-deficient"

  # solvable dimension
  inferdim = getFactorSolvableDim(dfg, fct, varid)
  # zdim = getFactorDim(fct)
  # vdim = getVariableDim(DFG.getVariable(dfg, varid))

  # TODO -- better to upsample before the projection
  Ndim = size(ptsbw,1)
  Npoints = size(ptsbw,2)
  # Assume we only have large particle population sizes, thanks to addNode!
  manis = getManifolds(dfg, varid)
  # manis = getSofttype(DFG.getVariable(dfg, varid)).manifolds # older
  p = AMP.manikde!(ptsbw, manis)
  if Npoints != N # this is where we control the overall particle set size
      p = resample(p,N)
  end
  return (p, inferdim)
end



"""
    $(SIGNATURES)

Add all potentials associated with this variable in `cliq` to `dens`.
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
    $SIGNATURES

Previous generation function converting LikelihoodMessage into Vector of densities before inference product.

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
  
  inferdim = 0.0
  # consolidate NBPMessages and potentials
  dens = Array{BallTreeDensity,1}()
  partials = Dict{Int, Vector{BallTreeDensity}}()
  wfac = Vector{Symbol}()
  # figure out which densities to use
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



## =============================================================================================
# Iterate over variables in clique
## =============================================================================================


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
    if true
      potprod = nothing
      densPts, inferdim = predictbelief(fgl, vsym, :, N=N, dbg=dbg, logger=logger)
    else
      # we'd like to do this more pre-emptive and then just execute -- just point and skip up only msgs
      densPts, potprod, inferdim = cliqGibbs(fgl, cliq, vsym, fmsgs, N, dbg, getSofttype(vert) |> getManifolds, logger)
    end

    if size(densPts,1)>0
      updvert = DFG.getVariable(fgl, vsym)  # TODO --  can we remove this duplicate getVert?
      setValKDE!(updvert, densPts, true, inferdim)
      # TODO perhpas more debugging inside `predictbelief`?
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
"""
function upGibbsCliqueDensity(dfg::AbstractDFG, cliq::TreeClique, 
                              inmsgs,
                              N::Int=100,
                              dbg::Bool=false,
                              iters::Int=3,
                              logger=ConsoleLogger()  ) where {T, T2}
  #
  with_logger(logger) do
    @info "up w $(length(inmsgs)) msgs"
  end
  # TODO -- some weirdness with: d,. = d = ., nothing
  mcmcdbg = Array{CliqGibbsMC,1}()
  d = Dict{Symbol,TreeBelief}()

  priorprods = Vector{CliqGibbsMC}()

  cliqdata = getCliqueData(cliq)


  # use nested structure for more efficient Chapman-Kolmogorov solution approximation
  if false
    IDS = [cliqdata.frontalIDs;cliqdata.separatorIDs] #inp.cliq.attributes["frontalIDs"]
    mcmcdbg, d = fmcmc!(dfg, cliq, inmsgs, IDS, N, iters, dbg, logger)
  else
    # NOTE -- previous mistake, must iterate over directsvarIDs also (or incorporate once at the right time)
    dummy, d = fmcmc!(dfg, cliq, inmsgs, cliqdata.directFrtlMsgIDs, N, 1, dbg, logger, true)
    if length(cliqdata.msgskipIDs) > 0
      dummy, dd = fmcmc!(dfg, cliq, inmsgs, cliqdata.msgskipIDs, N, 1, dbg, logger, true)
      for md in dd d[md[1]] = md[2]; end
    end
    if length(cliqdata.itervarIDs) > 0
      mcmcdbg, ddd = fmcmc!(dfg, cliq, inmsgs, cliqdata.itervarIDs, N, iters, dbg, logger, false)
      for md in ddd d[md[1]] = md[2]; end
    end
    if length(cliqdata.directPriorMsgIDs) > 0
      doids = setdiff(cliqdata.directPriorMsgIDs, cliqdata.msgskipIDs)
      priorprods, dddd = fmcmc!(dfg, cliq, inmsgs, doids, N, 1, dbg, logger, true)
      for md in dddd d[md[1]] = md[2]; end
    end
  end

  return d
end


## =============================================================================================
# Down solve
## =============================================================================================


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





#