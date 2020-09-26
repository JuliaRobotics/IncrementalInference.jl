
# starting to add exports here
export fetchCliqHistoryAll!


#global pidx
global pidx = 1
global pidl = 1
global pidA = 1
global thxl = nprocs() > 4 ? floor(Int,nprocs()*0.333333333) : 1

global WORKERPOOL = WorkerPool()

# upploc to control processes done local to this machine and separated from other
# highly loaded processes. upploc() should be used for dispatching short and burst
# of high bottle neck computations. Use upp2() for general multiple dispatch.
function upploc()
    global pidl, thxl
    N = nprocs()
    pidl = (thxl > 1 ? ((pidl < thxl && pidl != 0 && thxl > 2) ? (pidl+1)%(thxl+1) : 2) : (N == 1 ? 1 : 2))
    return pidl
end

# upp2() may refer to processes on a different machine. Do not associate upploc()
# with processes on a separate computer -- this will be become more complicated when
# remote processes desire their own short fast 'local' burst computations. We
# don't want all the data traveling back and forth for shorter jobs
function upp2()
  global pidx, thxl
  N = nprocs()
  pidx = (N > 1 ? ((pidx < N && pidx != 0) ? (pidx+1)%(N+1) : thxl+1) : 1) #2 -- thxl+1
  return pidx
end

function uppA()
  global pidA
  N = nprocs()
  pidA = (N > 1 ? ((pidA < N && pidA != 0) ? (pidA+1)%(N+1) : 2) : 1) #2 -- thxl+1
  return pidA
end

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
                                  dbg::Bool=false )::Float64
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
  return inferdim
end


function packFromLocalPartials!(fgl::G,
                                partials::Dict{Int, Vector{BallTreeDensity}},
                                cliq::TreeClique,
                                vsym::Symbol,
                                N::Int,
                                dbg::Bool=false  ) where G <: AbstractDFG
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

Multiply different dimensions from partial constraints individually.
"""
function productpartials!(pGM::Array{Float64,2},
                          dummy::BallTreeDensity,
                          partials::Dict{Int, Vector{BallTreeDensity}},
                          manis::T  )::Nothing where {T <: Tuple}
  #
  addopT, diffopT, getManiMu, getManiLam = buildHybridManifoldCallbacks(manis)
  for (dimnum,pp) in partials
    dummy2 = marginal(dummy,[dimnum])
    if length(pp) > 1
      pGM[dimnum,:], = prodAppxMSGibbsS(dummy2, pp, nothing, nothing, Niter=5, addop=(addopT[dimnum],), diffop=(diffopT[dimnum],), getMu=(getManiMu[dimnum],))
    else
      pGM[dimnum,:] = getPoints(pp[1])
    end
  end
  nothing
end

"""
    $(SIGNATURES)

Multiply various full and partial dimension constraints.
"""
function prodmultiplefullpartials( dens::Vector{BallTreeDensity},
                                   partials::Dict{Int, Vector{BallTreeDensity}},
                                   Ndims::Int,
                                   N::Int,
                                   manis::T ) where {T <: Tuple}
  #
  # TODO -- reuse memory rather than rand here
  pq = AMP.manifoldProduct(dens, manis, Niter=5)
  # pGM, = prodAppxMSGibbsS(dummy, dens, nothing, nothing, 5)
  for (dimnum,pp) in partials
    push!(pp, marginal(pq, [dimnum] ) )
    # push!(pp, AMP.manikde!(pGM[dimnum,:] ))
  end
  dummy = AMP.manikde!(rand(Ndims,N),[1.0], manis);
  pGM = getPoints(pq)
  productpartials!(pGM, dummy, partials, manis)
  return pGM
end

"""
    $(SIGNATURES)

Multiply a single full and several partial dimension constraints.
"""
function prodmultipleonefullpartials( dens::Vector{BallTreeDensity},
                                      partials::Dict{Int, Vector{BallTreeDensity}},
                                      Ndims::Int,
                                      N::Int,
                                      manis::T  ) where {T <: Tuple}
  #
  # TODO -- reuse memory rather than rand here
  # TODO -- should this be [1.0] or ones(Ndims)
  dummy = AMP.manikde!(rand(Ndims,N), [1.0], manis)
  denspts = getPoints(dens[1])
  pGM = deepcopy(denspts)
  for (dimnum,pp) in partials
    # @show (manis[dimnum],)
    push!(pp, AMP.manikde!(pGM[dimnum:dimnum,:], (manis[dimnum],) ))
  end
  productpartials!(pGM, dummy, partials, manis)
  return pGM
end

"""
    $SIGNATURES

Take product of `dens` accompanied by optional `partials` proposal belief densities.

Notes
-----
- `d` dimensional product approximation
- `partials` are treated as one dimensional

Future
------
- Incorporate ApproxManifoldProducts to process variables in individual batches.
"""
function productbelief(dfg::G,
                       vertlabel::Symbol,
                       dens::Vector{BallTreeDensity},
                       partials::Dict{Int, Vector{BallTreeDensity}},
                       N::Int;
                       dbg::Bool=false,
                       logger=ConsoleLogger()  ) where {G <: AbstractDFG}
  #
  vert = DFG.getVariable(dfg, vertlabel)
  manis = getSofttype(vert) |> getManifolds
  pGM = Array{Float64,2}(undef, 0,0)
  lennonp, lenpart = length(dens), length(partials)
  if lennonp > 1
    # multiple non-partials
    Ndims = Ndim(dens[1])
    with_logger(logger) do
      @info "[$(lennonp)x$(lenpart)p,d$(Ndims),N$(N)],"
    end
    pGM = prodmultiplefullpartials(dens, partials, Ndims, N, manis)
  elseif lennonp == 1 && lenpart >= 1
    # one non partial and one or more partials
    Ndims = Ndim(dens[1])
    with_logger(logger) do
      @info "[$(lennonp)x$(lenpart)p,d$(Ndims),N$(N)],"
    end
    pGM = prodmultipleonefullpartials(dens, partials, Ndims, N, manis)
  elseif lennonp == 0 && lenpart >= 1
    # only partials
    denspts = getPoints(getKDE(dfg, vertlabel))
    Ndims = size(denspts,1)
    with_logger(logger) do
      @info "[$(lennonp)x$(lenpart)p,d$(Ndims),N$(N)],"
    end
    dummy = AMP.manikde!(rand(Ndims,N), ones(Ndims), manis) # [1.0] # TODO -- reuse memory rather than rand here
    pGM = deepcopy(denspts)
    productpartials!(pGM, dummy, partials, manis)
  # elseif lennonp == 0 && lenpart == 1
  #   info("[prtl]")
  #   pGM = deepcopy(getVal(fg,vertid,api=localapi) )
  #   for (dimnum,pp) in partials
  #     pGM[dimnum,:] = getPoints(pp)
  #   end
  elseif lennonp == 1 && lenpart == 0
    # @info "[drct]"
    pGM = getPoints(dens[1])
  else
    with_logger(logger) do
      @warn "Unknown density product on variable=$(vert.label), lennonp=$(lennonp), lenpart=$(lenpart)"
    end
    pGM = Array{Float64,2}(undef, 0,1)
  end

  return pGM
end

"""
    $SIGNATURES

Compute the proposals of a destination vertex for each of `factors` and place the result
as belief estimates in both `dens` and `partials` respectively.

Notes
- TODO: also return if proposals were "dimension-deficient" (aka ~rank-deficient).
"""
function proposalbeliefs!(dfg::G,
                          destvertlabel::Symbol,
                          factors::Vector{F},
                          # inferddimproposal::Vector{Float64},
                          dens::Vector{BallTreeDensity},
                          partials::Dict{Int, Vector{BallTreeDensity}};
                          solveKey::Symbol=:default,
                          N::Int=100,
                          dbg::Bool=false)::Vector{Float64} where {G <: AbstractDFG, F <: DFGFactor}
  #
  inferddimproposal = Vector{Float64}()
  count = 0
  for fct in factors
    count += 1
    data = getSolverData(fct)
    p, inferd = findRelatedFromPotential(dfg, fct, destvertlabel, N, dbg, solveKey=solveKey)
    if data.fnc.partial   # partial density
      pardims = data.fnc.usrfnc!.partial
      for dimnum in pardims
        if haskey(partials, dimnum)
          push!(partials[dimnum], marginal(p,[dimnum]))
        else
          partials[dimnum] = BallTreeDensity[marginal(p,[dimnum])]
        end
      end
    else # full density
      push!(dens, p)
    end
    push!(inferddimproposal, inferd)
  end
  inferddimproposal
end

"""
    $SIGNATURES

Calculate the proposals and products on `destvert` using `factors` in factor graph `dfg`.

Notes
- Returns tuple of product and whether full dimensional (=true) or partial (=false).
"""
function predictbelief(dfg::AbstractDFG,
                       destvert::DFGVariable,
                       factors::Vector{<:DFGFactor};
                       N::Int=0,
                       dbg::Bool=false,
                       logger=ConsoleLogger()  ) # where {G <: , F <: }
  #
  destvertlabel = destvert.label
  dens = Array{BallTreeDensity,1}()
  partials = Dict{Int, Vector{BallTreeDensity}}()

  # determine number of particles to draw from the marginal
  nn = N != 0 ? N : size(getVal(destvert),2)

  # memory for if proposals are full dimension
  # inferdim = Float64[0.0;]

  # get proposal beliefs
  inferdim = proposalbeliefs!(dfg, destvertlabel, factors, dens, partials, N=nn, dbg=dbg)

  # take the product
  pGM = productbelief(dfg, destvertlabel, dens, partials, nn, dbg=dbg, logger=logger )

  return pGM, sum(inferdim)
end

function predictbelief(dfg::G,
                       destvertsym::Symbol,
                       factorsyms::Vector{Symbol};
                       N::Int=0,
                       dbg::Bool=false,
                       logger=ConsoleLogger()  ) where G <: AbstractDFG
  #
  factors = map(fsym -> DFG.getFactor(dfg, fsym), factorsyms)

  vert = DFG.getVariable(dfg, destvertsym)

  # determine the number of particles to draw from the marginal
  nn = N != 0 ? N : size(getVal(vert),2)

  # do the belief prediction
  predictbelief(dfg, vert, factors, N=nn, dbg=dbg, logger=logger)
end

function predictbelief(dfg::G,
                       destvertsym::Symbol,
                       factorsyms::Colon;
                       N::Int=0,
                       dbg::Bool=false,
                       logger=ConsoleLogger()) where G <: AbstractDFG
  #
  predictbelief(dfg, destvertsym, getNeighbors(dfg, destvertsym), N=N, dbg=dbg, logger=logger )
end

"""
    $(SIGNATURES)

Using factor graph object `dfg`, project belief through connected factors
(convolution with conditional) to variable `sym` followed by a approximate functional product.

Return: product belief, full proposals, partial dimension proposals, labels
"""
function localProduct(dfg::AbstractDFG,
                      sym::Symbol;
                      solveKey::Symbol=:default,
                      N::Int=100,
                      dbg::Bool=false,
                      logger=ConsoleLogger() )
  #
  # TODO -- converge this function with predictbelief for this node
  dens = Array{BallTreeDensity,1}()
  partials = Dict{Int, Vector{BallTreeDensity}}()
  lb = Symbol[]
  fcts = Vector{DFGFactor}()
  # vector of all neighbors as Symbols
  cf = getNeighbors(dfg, sym)
  for f in cf
    gfct = DFG.getFactor(dfg, f)
    push!(fcts, gfct)
    push!(lb, gfct.label)
  end

  # memory for if proposals are full dimension
  # inferdim = Vector{Float64}(undef, length(fcts))

  # get proposal beliefs
  inferdim = proposalbeliefs!(dfg, sym, fcts, dens, partials, solveKey=solveKey, N=N, dbg=dbg)

  # take the product
  pGM = productbelief(dfg, sym, dens, partials, N, dbg=dbg, logger=logger )
  vari = DFG.getVariable(dfg, sym)
  pp = AMP.manikde!(pGM, getSofttype(vari) |> getManifolds )

  return pp, dens, partials, lb, sum(inferdim)
end
localProduct(dfg::G, lbl::T; solveKey::Symbol=:default, N::Int=100, dbg::Bool=false) where {G <: AbstractDFG, T <: AbstractString} = localProduct(dfg, Symbol(lbl), solveKey=solveKey, N=N, dbg=dbg)


"""
    $(SIGNATURES)

Initialize the belief of a variable node in the factor graph struct.
"""
function initVariable!(fgl::G,
                       sym::Symbol;
                       N::Int=100 ) where G <: AbstractDFG
  #
  @warn "initVariable! has been displaced by doautoinit! or initManual! -- might be revived in the future"

  vert = getVariable(fgl, sym)
  belief,b,c,d,infdim  = localProduct(fgl, sym, N=N)
  setValKDE!(vert, belief)

  nothing
end


"""
    $(SIGNATURES)

Perform one step of the minibatch clique Gibbs operation for solving the Chapman-Kolmogov
trasit integral -- here involving separate approximate functional convolution and
product operations.
"""
function cliqGibbs(fg::G,
                   cliq::TreeClique,
                   vsym::Symbol,
                   inmsgs::Array{LikelihoodMessage,1},
                   N::Int,
                   dbg::Bool,
                   manis::T,
                   logger=ConsoleLogger()  ) where {G <: AbstractDFG, T <: Tuple}
  #
  # several optimizations can be performed in this function TODO

  # consolidate NBPMessages and potentials
  dens = Array{BallTreeDensity,1}()
  partials = Dict{Int, Vector{BallTreeDensity}}()
  wfac = Vector{Symbol}()

  inferdim = 0.0
  inferdim += packFromIncomingDensities!(dens, wfac, vsym, inmsgs, manis)
  inferdim += packFromLocalPotentials!(fg, dens, wfac, cliq, vsym, N)
  packFromLocalPartials!(fg, partials, cliq, vsym, N, dbg)

  potprod = !dbg ? nothing : PotProd(vsym, getVal(fg,vsym), Array{Float64,2}(undef, 0,0), dens, wfac)
      # pts,inferdim = predictbelief(dfg, vsym, useinitfct)  # for reference only
  pGM = productbelief(fg, vsym, dens, partials, N, dbg=dbg, logger=logger )
  if dbg  potprod.product = pGM  end

  # with_logger(logger) do
  #   @info "cliqGibbs -- end $vsym, inferdim=$inferdim, ls(fg)=$(ls(fg, vsym))"
  # end

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
function fmcmc!(fgl::G,
                cliq::TreeClique,
                fmsgs::Vector{LikelihoodMessage},
                lbls::Vector{Symbol},
                N::Int,
                MCMCIter::Int,
                dbg::Bool=false,
                logger=ConsoleLogger(),
                multithreaded::Bool=false  ) where G <: AbstractDFG
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

Calculate a fresh (single step) approximation to the variable `sym` in clique `cliq` as though during the upward message passing.  The full inference algorithm may repeatedly calculate successive apprimxations to the variables based on the structure of the clique, factors, and incoming messages.
Which clique to be used is defined by frontal variable symbols (`cliq` in this case) -- see `getClique(...)` for more details.  The `sym` symbol indicates which symbol of this clique to be calculated.  **Note** that the `sym` variable must appear in the clique where `cliq` is a frontal variable.
"""
function treeProductUp(fg::AbstractDFG,
                       tree::AbstractBayesTree,
                       cliq::Symbol,
                       sym::Symbol;
                       N::Int=100,
                       dbg::Bool=false  )
  #
  cliq = getClique(tree, cliq)
  cliqdata = getCliqueData(cliq)

  # get all the incoming (upward) messages from the tree cliques
  # convert incoming messages to Int indexed format (semi-legacy format)
  # upmsgssym = LikelihoodMessage[]
  # for cl in childCliqs(tree, cliq)
  #   msgdict = getUpMsgs(cl)
  #   dict = Dict{Symbol, TreeBelief}()
  #   for (dsy, btd) in msgdict.belief
  #     vari = getVariable(fg, dsy)
  #     # manis = getSofttype(vari).manifolds
  #     dict[dsy] = TreeBelief(btd.val, btd.bw, btd.inferdim, getSofttype(vari))
  #   end
  #   push!( upmsgssym, LikelihoodMessage(beliefDict=dict) )
  # end
  upmsgssym = getMsgsUpChildren(tree, cliq, TreeBelief)

  # perform the actual computation
  manis = getSofttype(getVariable(fg, sym)) |> getManifolds
  pGM, potprod, fulldim = cliqGibbs( fg, cliq, sym, upmsgssym, N, dbg, manis )

  return pGM, potprod
end


"""
    $(SIGNATURES)

Calculate a fresh---single step---approximation to the variable `sym` in clique `cliq` as though during the downward message passing.  The full inference algorithm may repeatedly calculate successive apprimxations to the variable based on the structure of variables, factors, and incoming messages to this clique.
Which clique to be used is defined by frontal variable symbols (`cliq` in this case) -- see `getClique(...)` for more details.  The `sym` symbol indicates which symbol of this clique to be calculated.  **Note** that the `sym` variable must appear in the clique where `cliq` is a frontal variable.
"""
function treeProductDwn(fg::G,
                        tree::AbstractBayesTree,
                        cliq::Symbol,
                        sym::Symbol;
                        N::Int=100,
                        dbg::Bool=false  ) where G <: AbstractDFG
  #
  @warn "treeProductDwn might not be working properly at this time. (post DFG v0.6 upgrade maintenance required)"
  cliq = getClique(tree, cliq)
  cliqdata = getCliqueData(cliq)

  # get the local variable id::Int identifier
  # vertid = fg.labelDict[sym]

  # get all the incoming (upward) messages from the tree cliques
  # convert incoming messages to Int indexed format (semi-legacy format)
  cl = parentCliq(tree, cliq)
  msgdict = getDwnMsgs(cl[1])
  dict = Dict{Int, TreeBelief}()
  for (dsy, btd) in msgdict
      dict[fg.IDs[dsy]] = TreeBelief(btd.val, btd.bw, btd.inferdim, getSofttype(getVariable(fg,sym)) )
  end
  dwnmsgssym = LikelihoodMessage[LikelihoodMessage(dict);]

  # perform the actual computation
  pGM, potprod, fulldim = cliqGibbs( fg, cliq, sym, dwnmsgssym, N, dbg ) #vertid

  return pGM, potprod, sym, dwnmsgssym
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



function dwnPrepOutMsg(fg::AbstractDFG,
                       cliq::TreeClique,
                       dwnMsgs::Array{LikelihoodMessage,1},
                       d::Dict{Symbol, T},
                       logger=ConsoleLogger()) where T
  # pack all downcoming conditionals in a dictionary too.
  with_logger(logger) do
    if cliq.index != 1 #TODO there may be more than one root
      @info "Dwn msg keys $(keys(dwnMsgs[1].belief))"
      @info "fg vars $(ls(fg))"
    end # ignore root, now incoming dwn msg
  end
  m = LikelihoodMessage()
  i = 0
  for vid in getCliqueData(cliq).frontalIDs
    m.belief[vid] = deepcopy(d[vid]) # TODO -- not sure if deepcopy is required
  end
  for cvid in getCliqueData(cliq).separatorIDs
    i+=1
    # TODO -- convert to points only since kde replace by rkhs in future
    m.belief[cvid] = deepcopy(dwnMsgs[1].belief[cvid]) # TODO -- maybe this can just be a union(,)
  end
  return m
end

"""
    $SIGNATURES

Perform Chapman-Kolmogorov transit integral approximation for `cliq` in downward pass direction.

Notes
- Only update frontal variables of the clique.
"""
function downGibbsCliqueDensity(fg::G,
                                cliq::TreeClique,
                                dwnMsgs::Array{LikelihoodMessage,1},
                                N::Int=100,
                                MCMCIter::Int=3,
                                dbg::Bool=false,
                                usemsgpriors::Bool=false,
                                logger=ConsoleLogger()) where G <: AbstractDFG
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
function downGibbsCliqueDensity(fg::G,
                                cliq::TreeClique,
                                dwnMsgs::LikelihoodMessage,
                                N::Int=100,
                                MCMCIter::Int=3,
                                dbg::Bool=false,
                                usemsgpriors::Bool=false,
                                logger=ConsoleLogger()) where G <: AbstractDFG
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

Set the color of a cliq in the Bayes (Junction) tree.
"""
function setCliqDrawColor(cliq::TreeClique, fillcolor::String)::Nothing
  cliq.attributes["fillcolor"] = fillcolor
  cliq.attributes["style"] = "filled"
  nothing
end

function getCliqueDrawColor(cliq::TreeClique)
  haskey(cliq.attributes, "fillcolor") ? cliq.attributes["fillcolor"] : nothing
end

"""
    $(SIGNATURES)

Update cliq `cliqID` in Bayes (Juction) tree `bt` according to contents of `ddt` -- intended use is to update main clique after a downward belief propagation computation has been completed per clique.
"""
function updateFGBT!(fg::AbstractDFG,
                     bt::AbstractBayesTree,
                     cliqID::Int,
                     drt::DownReturnBPType;
                     dbg::Bool=false,
                     fillcolor::String="",
                     logger=ConsoleLogger()  )
    #
    cliq = getClique(bt, cliqID)
    # if dbg
    #   cliq.attributes["debugDwn"] = deepcopy(drt.dbgDwn)
    # end
    setDwnMsg!(cliq, drt.keepdwnmsgs)
    # TODO move to drawTree
    if fillcolor != ""
      setCliqDrawColor(cliq, fillcolor)
    end
    with_logger(logger) do
      for (sym, emsg) in drt.IDvals
        #TODO -- should become an update call
        updvert = DFG.getVariable(fg, sym)
        # TODO -- not sure if deepcopy is required , updatePPE=true)
        @info "updateFGBT, DownReturnBPType, sym=$sym, current inferdim val=$(getVariableInferredDim(updvert))"
        setValKDE!(updvert, deepcopy(emsg) )
        updvert = DFG.getVariable(fg, sym)
        @info "updateFGBT, DownReturnBPType, sym=$sym, inferdim=$(emsg.inferdim), newval=$(getVariableInferredDim(updvert))"
      end
    end
    nothing
end

"""
    $(SIGNATURES)

Update cliq `cliqID` in Bayes (Juction) tree `bt` according to contents of `urt` -- intended use is to update main clique after a upward belief propagation computation has been completed per clique.
"""
function updateFGBT!(fg::AbstractDFG,
                     cliq::TreeClique,
                     IDvals::Dict{Symbol, TreeBelief};
                     dbg::Bool=false,
                     fillcolor::String="",
                     logger=ConsoleLogger()  )
  #
  # if dbg
  #   # TODO find better location for the debug information (this is old code)
  #   cliq.attributes["debug"] = deepcopy(urt.dbgUp)
  # end
  if fillcolor != ""
    setCliqDrawColor(cliq, fillcolor)
  end
  for (id,dat) in IDvals
    with_logger(logger) do
      @info "updateFGBT! up -- update $id, inferdim=$(dat.inferdim)"
    end
    updvert = DFG.getVariable(fg, id)
    setValKDE!(updvert, deepcopy(dat), true) ## TODO -- not sure if deepcopy is required
  end
  with_logger(logger) do
    @info "updateFGBT! up -- updated $(getLabel(cliq))"
  end
  nothing
end




"""
    $SIGNATURES

Approximate Chapman-Kolmogorov transit integral and return separator marginals as messages to pass up the Bayes (Junction) tree, along with additional clique operation values for debugging.

Notes
- `onduplicate=true` by default internally uses deepcopy of factor graph and Bayes tree, and does **not** update the given objects.  Set false to update `fgl` and `treel` during compute.

Future
- TODO: internal function chain is too long and needs to be refactored for maintainability.
"""
function approxCliqMarginalUp!(csmc::CliqStateMachineContainer,
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

Special internal function to try return the clique data if succesfully identified in `othertree::AbstractBayesTree`,
based on contents of `seeksSimilar::BayesTreeNodeData`.

Notes
- Used to identify and skip similar cliques (i.e. recycle computations)
"""
function attemptTreeSimilarClique(othertree::AbstractBayesTree, seeksSimilar::BayesTreeNodeData)::TreeClique
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
  return otherCliq
end



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
    error(err)
  end
  # if !(clst in [:upsolved; :downsolved; :marginalized])
  #   error("Clique $(cliq.index), initInferTreeUp! -- cliqInitSolveUp! did not arrive at the desired solution statu: $clst")
  # end
  return history
end


"""
    $SIGNATURES

Fetch solver history from clique state machines that have completed their async Tasks and store in the `hist::Dict{Int,Tuple}` dictionary.
"""
function fetchCliqHistoryAll!(smt::Vector{Task},
                              hist::Dict{Int,Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}=Dict{Int,Vector{Tuple{DateTime, Int,
                                                                    Function, CliqStateMachineContainer}}}() )
  #
  for i in 1:length(smt)
    sm = smt[i]
    # only fetch states that have completed processing
    if sm.state == :done
      haskey(hist, i) ? @warn("overwriting existing history key $i") : nothing
      hist[i] = fetch(sm)
    end
  end
  hist
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
function asyncTreeInferUp!(dfg::AbstractDFG,
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
