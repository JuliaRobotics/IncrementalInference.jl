#global pidx
global pidx = 1
global pidl = 1
global pidA = 1
global thxl = nprocs() > 4 ? floor(Int,nprocs()*0.333333333) : 1

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


function packFromIncomingDensities!(dens::Vector{BallTreeDensity},
                                    wfac::Vector{Symbol},
                                    vsym::Symbol,
                                    inmsgs::Array{NBPMessage,1},
                                    manis::T )::Float64 where {T <: Tuple}
  #
  inferdim = 0.0
  for m in inmsgs
    for psym in keys(m.p)
      if psym == vsym
        pdi = m.p[vsym] # ::EasyMessage
        push!(dens, manikde!(pdi.pts, pdi.bws, pdi.manifolds) ) # kde!(pdi.pts, pdi.bws)
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
function packFromLocalPotentials!(dfg::G,
                                  dens::Vector{BallTreeDensity},
                                  wfac::Vector{Symbol},
                                  cliq::Graphs.ExVertex,
                                  vsym::Symbol,
                                  N::Int,
                                  dbg::Bool=false )::Float64 where G <: AbstractDFG
  #
  inferdim = 0.0
  for idfct in getData(cliq).potentials
    fct = DFG.getFactor(dfg, idfct)
    data = getData(fct)
    # skip partials here, will be caught in packFromLocalPartials!
    if length( findall(data.fncargvID .== vsym) ) >= 1 && !data.fnc.partial
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
                                cliq::Graphs.ExVertex,
                                vsym::Symbol,
                                N::Int,
                                dbg::Bool=false  ) where G <: AbstractDFG
  #

  for idfct in getData(cliq).potentials
    vert = DFG.getFactor(fgl, idfct)
    data = getData(vert)
    if length( findall(data.fncargvID .== vsym) ) >= 1 && data.fnc.partial
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
  manis = getSofttype(vert).manifolds
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
                          N::Int=100,
                          dbg::Bool=false)::Vector{Float64} where {G <: AbstractDFG, F <: DFGFactor}
  #
  inferddimproposal = Vector{Float64}()
  count = 0
  for fct in factors
    count += 1
    data = getData(fct)
    p, inferd = findRelatedFromPotential(dfg, fct, destvertlabel, N, dbg)
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
function predictbelief(dfg::G,
                       destvert::DFGVariable,
                       factors::Vector{F};
                       N::Int=0,
                       dbg::Bool=false,
                       logger=ConsoleLogger()  ) where {G <: AbstractDFG, F <: DFGNode}
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
function localProduct(dfg::G,
                      sym::Symbol;
                      N::Int=100,
                      dbg::Bool=false,
                      logger=ConsoleLogger() ) where G <: AbstractDFG
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
  inferdim = proposalbeliefs!(dfg, sym, fcts, dens, partials, N=N, dbg=dbg)

  # take the product
  pGM = productbelief(dfg, sym, dens, partials, N, dbg=dbg, logger=logger )
  vari = DFG.getVariable(dfg, sym)
  pp = AMP.manikde!(pGM, getSofttype(vari).manifolds )

  return pp, dens, partials, lb, sum(inferdim)
end
localProduct(dfg::G, lbl::T; N::Int=100, dbg::Bool=false) where {G <: AbstractDFG, T <: AbstractString} = localProduct(dfg, Symbol(lbl), N=N, dbg=dbg)


"""
    $(SIGNATURES)

Initialize the belief of a variable node in the factor graph struct.
"""
function initVariable!(fgl::G,
                       sym::Symbol;
                       N::Int=100 ) where G <: AbstractDFG
  #
  @warn "initVariable! has been displaced by doautoinit! or manualinit! -- might be revived in the future"

  vert = getVariable(fgl, sym)
  belief,b,c,d,infdim  = localProduct(fgl, sym)
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
                   cliq::Graphs.ExVertex,
                   vsym::Symbol,
                   inmsgs::Array{NBPMessage,1},
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
function compileFMCMessages(fgl::G, lbls::Vector{Symbol}, logger=ConsoleLogger()) where G <: AbstractDFG
  d = Dict{Symbol,EasyMessage}()
  for vsym in lbls
    vert = DFG.getVariable(fgl,vsym)
    pden = getKDE(vert)
    bws = vec(getBW(pden)[:,1])
    manis = getSofttype(vert).manifolds
    d[vsym] = EasyMessage(getVal(vert), bws, manis, getData(vert).inferdim)
    with_logger(logger) do
      @info "fmcmc! -- getData(vert=$(vert.label)).inferdim=$(getData(vert).inferdim)"
    end
  end
  return d
end

function doFMCIteration(fgl::G,
                        vsym::Symbol,
                        cliq::Graphs.ExVertex,
                        fmsgs,
                        N::Int,
                        dbg::Bool,
                        logger=ConsoleLogger()  ) where G <: AbstractDFG
  #
  vert = DFG.getVariable(fgl, vsym)
  if !getData(vert).ismargin
    # we'd like to do this more pre-emptive and then just execute -- just point and skip up only msgs
    densPts, potprod, inferdim = cliqGibbs(fgl, cliq, vsym, fmsgs, N, dbg, getSofttype(vert).manifolds, logger)

    if size(densPts,1)>0
      updvert = DFG.getVariable(fgl, vsym)  # TODO --  can we remove this duplicate getVert?
      setValKDE!(updvert, densPts, true, inferdim)
      if dbg
        push!(dbgvals.prods, potprod)
        push!(dbgvals.lbls, Symbol(updvert.label))
      end
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
                cliq::Graphs.ExVertex,
                fmsgs::Vector{NBPMessage},
                lbls::Vector{Symbol},
                N::Int,
                MCMCIter::Int,
                dbg::Bool=false,
                logger=ConsoleLogger(),
                multithreaded::Bool=false  ) where G <: AbstractDFG
  #
  with_logger(logger) do
    @info "---------- successive fnc approx ------------$(cliq.attributes["label"])"
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

      # vert = DFG.getVariable(fgl, vsym)
      # if !getData(vert).ismargin
      #   # we'd like to do this more pre-emptive and then just execute -- just point and skip up only msgs
      #   densPts, potprod, inferdim = cliqGibbs(fgl, cliq, vsym, fmsgs, N, dbg, getSofttype(vert).manifolds, logger)
      #
      #   if size(densPts,1)>0
      #     updvert = DFG.getVariable(fgl, vsym)  # TODO --  can we remove this duplicate getVert?
      #     setValKDE!(updvert, densPts, true, inferdim)
      #     if dbg
      #       push!(dbgvals.prods, potprod)
      #       push!(dbgvals.lbls, Symbol(updvert.label))
      #     end
      #   end
      # end
    end
    !dbg ? nothing : push!(mcmcdbg, dbgvals)
  end

  # populate dictionary for return NBPMessage in multiple dispatch
  msgdict = compileFMCMessages(fgl, lbls, logger)

  return mcmcdbg, msgdict
end

"""
    $SIGNATURES

MUST BE REFACTORED OR DEPRECATED.  Seems like a wasteful function.
"""
function upPrepOutMsg!(d::Dict{Symbol,EasyMessage}, IDs::Vector{Symbol}) #Array{Float64,2}
  @info "Outgoing msg density on: "
  len = length(IDs)
  m = NBPMessage(Dict{Symbol,EasyMessage}())
  for id in IDs
    m.p[id] = d[id]
  end
  return m
end

"""
    $(SIGNATURES)

Calculate a fresh (single step) approximation to the variable `sym` in clique `cliq` as though during the upward message passing.  The full inference algorithm may repeatedly calculate successive apprimxations to the variables based on the structure of the clique, factors, and incoming messages.
Which clique to be used is defined by frontal variable symbols (`cliq` in this case) -- see `whichCliq(...)` for more details.  The `sym` symbol indicates which symbol of this clique to be calculated.  **Note** that the `sym` variable must appear in the clique where `cliq` is a frontal variable.
"""
function treeProductUp(fg::G,
                       tree::BayesTree,
                       cliq::Symbol,
                       sym::Symbol;
                       N::Int=100,
                       dbg::Bool=false  ) where G <: AbstractDFG
  #
  cliq = whichCliq(tree, cliq)
  cliqdata = getData(cliq)
  # IDS = [cliqdata.frontalIDs; cliqdata.conditIDs]

  # get the local variable id::Int identifier
  vertid = fg.IDs[sym]

  # get all the incoming (upward) messages from the tree cliques
  # convert incoming messages to Int indexed format (semi-legacy format)
  upmsgssym = NBPMessage[]
  for cl in childCliqs(tree, cliq)
    msgdict = upMsg(cl)
    dict = Dict{Int, EasyMessage}()
    for (dsy, btd) in msgdict
      manis = getSofttype(getVert(fg, dsy, api=localapi)).manifolds

      dict[fg.IDs[dsy]] = convert(EasyMessage, btd, manis)
    end
    push!( upmsgssym, NBPMessage(dict) )
  end

  # perform the actual computation
  manis = getSofttype(getVert(fg, vertid, api=localapi)).manifolds
  pGM, potprod, fulldim = cliqGibbs( fg, cliq, vertid, upmsgssym, N, dbg, manis )

  return pGM, potprod
end


"""
    $(SIGNATURES)

Calculate a fresh---single step---approximation to the variable `sym` in clique `cliq` as though during the downward message passing.  The full inference algorithm may repeatedly calculate successive apprimxations to the variable based on the structure of variables, factors, and incoming messages to this clique.
Which clique to be used is defined by frontal variable symbols (`cliq` in this case) -- see `whichCliq(...)` for more details.  The `sym` symbol indicates which symbol of this clique to be calculated.  **Note** that the `sym` variable must appear in the clique where `cliq` is a frontal variable.
"""
function treeProductDwn(fg::G,
                        tree::BayesTree,
                        cliq::Symbol,
                        sym::Symbol;
                        N::Int=100,
                        dbg::Bool=false  ) where G <: AbstractDFG
  #
  cliq = whichCliq(tree, cliq)
  cliqdata = getData(cliq)

  # get the local variable id::Int identifier
  vertid = fg.IDs[sym]

  # get all the incoming (upward) messages from the tree cliques
  # convert incoming messages to Int indexed format (semi-legacy format)
  cl = parentCliq(tree, cliq)
  msgdict = getDwnMsgs(cl[1])
  dict = Dict{Int, EasyMessage}()
  for (dsy, btd) in msgdict
      dict[fg.IDs[dsy]] = convert(EasyMessage, btd)
  end
  dwnmsgssym = NBPMessage[NBPMessage(dict);]

  # perform the actual computation
  pGM, potprod, fulldim = cliqGibbs( fg, cliq, vertid, dwnmsgssym, N, dbg )

  return pGM, potprod, vertid, dwnmsgssym
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
  d = Dict{Symbol,EasyMessage}()

  priorprods = Vector{CliqGibbsMC}()

  cliqdata = getData(inp.cliq)

  with_logger(logger) do
    for el in inp.sendmsgs, (id,msg) in el.p
      @info "inp.sendmsgs[$id].inferdim=$(msg.inferdim)"
    end
  end

  # use nested structure for more fficient Chapman-Kolmogorov solution approximation
  if false
    IDS = [cliqdata.frontalIDs;cliqdata.conditIDs] #inp.cliq.attributes["frontalIDs"]
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

  #m = upPrepOutMsg!(inp.fg, inp.cliq, inp.sendmsgs, condids, N)
  m = upPrepOutMsg!(d, cliqdata.conditIDs)

  # TODO can remove this outmsglbl Symbol => Symbol
  outmsglbl = Dict{Symbol, Symbol}()
  if dbg
    for (ke, va) in m.p
      outmsglbl[Symbol(inp.fg.g.vertices[ke].label)] = ke
    end
  end

  # prepare and convert upward belief messages
  upmsgs = TempBeliefMsg() #Dict{Symbol, BallTreeDensity}()
  # @show collect(keys(inp.fg.g.vertices))
  for (msgsym, val) in m.p
    upmsgs[msgsym] = convert(Tuple{BallTreeDensity,Float64}, val) # (convert(BallTreeDensity, val), getVariableInferredDim(inp.fg,msgsym))
  end
  setUpMsg!(inp.cliq, upmsgs)

  # flag cliq as definitely being initialized
  cliqdata.upsolved = true

  mdbg = !dbg ? DebugCliqMCMC() : DebugCliqMCMC(mcmcdbg, m, outmsglbl, priorprods)
  return UpReturnBPType(m, mdbg, d, upmsgs, true)
end



function dwnPrepOutMsg(fg::G,
                       cliq::Graphs.ExVertex,
                       dwnMsgs::Array{NBPMessage,1},
                       d::Dict{Symbol, T},
                       logger=ConsoleLogger()) where {G <: AbstractDFG, T}
  # pack all downcoming conditionals in a dictionary too.
  with_logger(logger) do
    if cliq.index != 1
      @info "Dwn msg keys $(keys(dwnMsgs[1].p))"
      @info "fg vars $(ls(fg))"
    end # ignore root, now incoming dwn msg
  end
  m = NBPMessage(Dict{Symbol,T}())
  i = 0
  for vid in getData(cliq).frontalIDs
    m.p[vid] = deepcopy(d[vid]) # TODO -- not sure if deepcopy is required
  end
  for cvid in getData(cliq).conditIDs
    i+=1
    # TODO -- convert to points only since kde replace by rkhs in future
    m.p[cvid] = deepcopy(dwnMsgs[1].p[cvid]) # TODO -- maybe this can just be a union(,)
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
                                cliq::Graphs.ExVertex,
                                dwnMsgs::Array{NBPMessage,1},
                                N::Int=100,
                                MCMCIter::Int=3,
                                dbg::Bool=false,
                                usemsgpriors::Bool=false,
                                logger=ConsoleLogger()) where G <: AbstractDFG
  #
  # TODO standardize function call to have similar stride to upGibbsCliqueDensity
  # @info "down"
  with_logger(logger) do
    @info "cliq=$(cliq.index), downGibbsCliqueDensity -- going for down fmcmc, keys=$(keys(dwnMsgs))"
  end
  fmmsgs = usemsgpriors ? Array{NBPMessage,1}() : dwnMsgs
  frtls = getFrontals(cliq)

  # TODO, do better check if there is structure between multiple frontals
  niters = length(frtls) == 1 ? 1 : MCMCIter
  # TODO standize with upsolve and variable solver order
  mcmcdbg, d = fmcmc!(fg, cliq, fmmsgs, frtls, N, niters, dbg)
  m = dwnPrepOutMsg(fg, cliq, dwnMsgs, d, logger)

  outmsglbl = Dict{Symbol, Int}()
  if dbg
    for (ke, va) in m.p
      outmsglbl[Symbol(fg.g.vertices[ke].label)] = ke
    end
  end

  with_logger(logger) do
    @info "cliq=$(cliq.index), downGibbsCliqueDensity -- convert to BallTreeDensities."
  end

  # Always keep dwn messages in cliq data
  dwnkeepmsgs = TempBeliefMsg() # Dict{Symbol, BallTreeDensity}()
  for (msgsym, val) in m.p
    dwnkeepmsgs[msgsym] = convert(Tuple{BallTreeDensity,Float64}, val)
  end
  setDwnMsg!(cliq, dwnkeepmsgs)

  # down solving complete, set flag
  getData(cliq).downsolved = true

  mdbg = !dbg ? DebugCliqMCMC() : DebugCliqMCMC(mcmcdbg, m, outmsglbl, CliqGibbsMC[])
  with_logger(logger) do
    @info "cliq=$(cliq.index), downGibbsCliqueDensity -- finished."
  end
  return DownReturnBPType(m, mdbg, d, dwnkeepmsgs)
end
function downGibbsCliqueDensity(fg::G,
                                cliq::Graphs.ExVertex,
                                dwnMsgs::TempBeliefMsg, # Dict{Symbol,BallTreeDensity},
                                N::Int=100,
                                MCMCIter::Int=3,
                                dbg::Bool=false,
                                usemsgpriors::Bool=false,
                                logger=ConsoleLogger()) where G <: AbstractDFG
  #
  with_logger(logger) do
    @info "cliq=$(cliq.index), downGibbsCliqueDensity -- convert BallTreeDensities to NBPMessages."
  end
  ind = Dict{Symbol, EasyMessage}()
  sflbls = getVariableIds(fg)
  for (lbl, bel) in dwnMsgs
	  if lbl in sflbls
	    ind[lbl] = convert(EasyMessage, bel, getManifolds(fg, lbl))
    end
  end
  ndms = NBPMessage[NBPMessage(ind);]
  with_logger(logger) do
    @info "cliq=$(cliq.index), downGibbsCliqueDensity -- call with NBPMessages."
  end
  downGibbsCliqueDensity(fg, cliq, ndms, N, MCMCIter, dbg, usemsgpriors, logger)
end

"""
    $SIGNATURES

Set the color of a cliq in the Bayes (Junction) tree.
"""
function setCliqDrawColor(cliq::Graphs.ExVertex, fillcolor::String)::Nothing
  cliq.attributes["fillcolor"] = fillcolor
  cliq.attributes["style"] = "filled"
  nothing
end

"""
    $(SIGNATURES)

Update cliq `cliqID` in Bayes (Juction) tree `bt` according to contents of `ddt` -- intended use is to update main clique after a downward belief propagation computation has been completed per clique.
"""
function updateFGBT!(fg::G,
                     bt::BayesTree,
                     cliqID::Int,
                     drt::DownReturnBPType;
                     dbg::Bool=false,
                     fillcolor::String="",
                     logger=ConsoleLogger()  ) where G <: AbstractDFG
    #
    cliq = bt.cliques[cliqID]
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
        # TODO -- not sure if deepcopy is required , updateMAPest=true)
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
function updateFGBT!(fg::G,
                     cliq::Graphs.ExVertex,
                     urt::UpReturnBPType;
                     dbg::Bool=false,
                     fillcolor::String="",
                     logger=ConsoleLogger()  ) where G <: AbstractDFG
  #
  if dbg
    cliq.attributes["debug"] = deepcopy(urt.dbgUp)
  end
  setUpMsg!(cliq, urt.keepupmsgs)
  # move to drawTree
  if fillcolor != ""
    setCliqDrawColor(cliq, fillcolor)
  end
  # cliqFulldim = true
  for (id,dat) in urt.IDvals
    # cliqFulldim &= dat.fulldim
    with_logger(logger) do
      @info "updateFGBT! up -- update $id, inferdim=$(dat.inferdim)"
    end
    updvert = DFG.getVariable(fg, id)
    setValKDE!(updvert, deepcopy(dat), true) ## TODO -- not sure if deepcopy is required
  end
  with_logger(logger) do
    @info "updateFGBT! up -- updated $(cliq.attributes["label"])"
  end
  nothing
end

function updateFGBT!(fg::G,
                     bt::BayesTree,
                     cliqID::Int,
                     urt::UpReturnBPType;
                     dbg::Bool=false, fillcolor::String=""  ) where G <: AbstractDFG
  #
  cliq = bt.cliques[cliqID]
  cliq = bt.cliques[cliqID]
  updateFGBT!( fg, cliq, urt, dbg=dbg, fillcolor=fillcolor )
end


"""
    $SIGNATURES

Get and return upward belief messages as stored in child cliques from `treel::BayesTree`.

Notes
- Use last parameter to select the return format.
"""
function getCliqChildMsgsUp(fg_::G,
                            treel::BayesTree,
                            cliq::Graphs.ExVertex,
                            ::Type{EasyMessage} ) where G <: AbstractDFG
  #
  childmsgs = NBPMessage[]
  for child in getChildren(treel, cliq)
    nbpchild = NBPMessage(Dict{Symbol,EasyMessage}())
    for (key, bel) in getUpMsgs(child)
      manis = getManifolds(fg_, key)
      # inferdim = getVariableInferredDim(fg_, key)
      nbpchild.p[key] = convert(EasyMessage, bel, manis)
    end
    push!(childmsgs, nbpchild)
  end
  return childmsgs
end

function getCliqChildMsgsUp(treel::BayesTree, cliq::Graphs.ExVertex, ::Type{BallTreeDensity})
  childmsgs = Dict{Symbol,Vector{Tuple{BallTreeDensity,Float64}}}()  # Vector{Bool}
  for child in getChildren(treel, cliq)
    for (key, bel) in getUpMsgs(child)
      # id = fg_.IDs[key]
      # manis = getManifolds(fg_, id)
      if !haskey(childmsgs, key)
        childmsgs[key] = Vector{Tuple{BallTreeDensity, Float64}}()  # Vector{Bool}
      end
      push!(childmsgs[key], bel )
    end
  end
  return childmsgs
end

"""
    $SIGNATURES

Get the latest down message from the parent node (without calculating anything).

Notes
- Different from down initialization messages that do calculate new values -- see `prepCliqInitMsgsDown!`.
- Basically converts function `getDwnMsgs` from `Dict{Symbol,BallTreeDensity}` to `Dict{Symbol,Vector{BallTreeDensity}}`.
"""
function getCliqParentMsgDown(treel::BayesTree, cliq::Graphs.ExVertex)
  downmsgs = Dict{Symbol,Vector{Tuple{BallTreeDensity, Float64}}}()
  for prnt in getParent(treel, cliq)
    for (key, bel) in getDwnMsgs(prnt)
      if !haskey(downmsgs, key)
        downmsgs[key] = Vector{Tuple{BallTreeDensity, Float64}}()
      end
      # TODO insert true inferred dim
      push!(downmsgs[key], bel)
    end
  end
  return downmsgs
end



"""
    $SIGNATURES

Approximate Chapman-Kolmogorov transit integral and return separator marginals as messages to pass up the Bayes (Junction) tree, along with additional clique operation values for debugging.

Notes
- `onduplicate=true` by default internally uses deepcopy of factor graph and Bayes tree, and does **not** update the given objects.  Set false to update `fgl` and `treel` during compute.

Future
- TODO: internal function chain is too long and needs to be refactored for maintainability.
"""
function approxCliqMarginalUp!(fgl::G,
                               treel::BayesTree,
                               csym::Symbol,
                               onduplicate=true;
                               N::Int=100,
                               dbg::Bool=false,
                               iters::Int=3,
                               drawpdf::Bool=false,
                               multiproc::Bool=true,
                               logger=ConsoleLogger()  ) where G <: AbstractDFG
  #
  fg_ = onduplicate ? deepcopy(fgl) : fgl
  onduplicate
  with_logger(logger) do
    @warn "rebuilding new Bayes tree on deepcopy of factor graph"
  end
  tree_ = onduplicate ? wipeBuildNewTree!(fgl) : treel


  # copy up and down msgs that may already exists
  if onduplicate
    for (id, cliq) in tree_.cliques
      setUpMsg!(tree_.cliques[cliq.index], getUpMsgs(cliq))
      setDwnMsg!(tree_.cliques[cliq.index], getDwnMsgs(cliq))
    end
  end

  cliq = whichCliq(tree_, csym)
  # setCliqDrawColor(cliq, "red")

  # get incoming cliq messaged upward from child cliques
  childmsgs = getCliqChildMsgsUp(fg_, tree_, cliq, EasyMessage)

  # TODO use subgraph copy of factor graph for operations and transfer frontal variables only

  with_logger(logger) do
    @info "=== start Clique $(cliq.attributes["label"]) ======================"
  end
  ett = FullExploreTreeType(fg_, nothing, cliq, nothing, childmsgs)
  urt = UpReturnBPType()
  if multiproc
    cliqc = deepcopy(cliq)
    cliqcd = getData(cliqc)
    # redirect to new unused so that CAN be serialized
    cliqcd.initUpChannel = Channel{Symbol}(1)
    cliqcd.initDownChannel = Channel{Symbol}(1)
    cliqcd.solveCondition = Condition()
    cliqcd.statehistory = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}()
    ett.cliq = cliqc
    # TODO create new dedicate file for separate process to log with
    urt = remotecall_fetch(upGibbsCliqueDensity, upp2(), ett, N, dbg, iters)
  else
    with_logger(logger) do
      @info "Single process upsolve clique=$(cliq.index)"
    end
    urt = upGibbsCliqueDensity(ett, N, dbg, iters, logger)
  end

  # is clique fully upsolved or only partially?
  with_logger(logger) do
    updateFGBT!(fgl, cliq, urt, dbg=dbg, fillcolor="brown", logger=logger)
  end

  # is clique fulldim/partially solved?
  # for (id, val) in urt.IDvals
  #     cliqFulldim += val.fulldim
  # end

  # set clique color accordingly, using local memory
  setCliqDrawColor(cliq, isCliqFullDim(fg_, cliq) ? "pink" : "tomato1")

  drawpdf ? drawTree(tree_) : nothing
  with_logger(logger) do
    @info "=== end Clique $(cliq.attributes["label"]) ========================"
  end
  return urt
end

"""
    $SIGNATURES

Approximate Chapman-Kolmogorov transit integral and return separator marginals as messages to pass up the Bayes (Junction) tree, along with additional clique operation values for debugging.

Notes
=====
- `onduplicate=true` by default internally uses deepcopy of factor graph and Bayes tree, and does **not** update the given objects.  Set false to update `fgl` and `treel` during compute.
"""
function doCliqInferenceUp!(fgl::FactorGraph,
                            treel::BayesTree,
                            csym::Symbol,
                            onduplicate=true;
                            N::Int=100,
                            dbg::Bool=false,
                            iters::Int=3,
                            drawpdf::Bool=false,
                            multiproc::Bool=true,
                            logger=ConsoleLogger()   )
  #
  approxCliqMarginalUp!(fgl, treel, csym, onduplicate; N=N, dbg=dbg, iters=iters, drawpdf=drawpdf, multiproc=multiproc, logger=logger  )
end

# """
#     $(TYPEDSIGNATURES)
#
# Perform Chapman-Kolmogorov transit integral procedure for a given clique, specifically for the upward direction.
#
# Notes:
# -----
# * Assumes the same procedure has completed for the child cliques, since upward messages are required from them.
# * Can adjust the number of `iters::Int=3` must be performed on the `itervars` of this clique.
# """
# function doCliqInferenceUp!(fgl::FactorGraph,
#                             treel::BayesTree,
#                             cliql::Graphs.ExVertex;
#                             N::Int=100,
#                             dbg::Bool=false,
#                             iters::Int=3  )
#   #
#   # get children
#   childr = childCliqs(treel, cliql)
#
#   # get upward messages from children
#   upmsgs = NBPMessage[]
#   for child in childr
#     frsym = getSym(fgl, getFrontals(child)[1])
#     ret = getUpMsgs(treel, frsym)
#     @show typeof(ret)
#     newmsg = NBPMessage()
#     for (id, val) in ret
#       newmsg[id] = convert(EasyMessage, val, manis)
#     end
#     push!(upmsgs, newmsg)
#   end
#
#   ett = ExploreTreeType(fgl, treel, cliql, nothing, upmsgs)
#
#   urt = upGibbsCliqueDensity(ett, N, dbg, iters)
#   return urt.keepupmsgs
# end




## NOTE REMOVED MANY RECURSIVE OR ITERATIVE STACK FUNCTIONS HERE



"""
    $SIGNATURES

Set all up `upsolved` and `downsolved` cliq data flags `to::Bool=false`.
"""
function setAllSolveFlags!(treel::BayesTree, to::Bool=false)::Nothing
  for (id, cliq) in treel.cliques
    cliqdata = getData(cliq)
    cliqdata.initialized = :null
    cliqdata.upsolved = to
    cliqdata.downsolved = to
  end
  nothing
end

"""
    $SIGNATURES

Return true or false depending on whether the tree has been fully initialized/solved/marginalized.
"""
function isTreeSolved(treel::BayesTree; skipinitialized::Bool=false)::Bool
  acclist = Symbol[:upsolved; :downsolved; :marginalized]
  skipinitialized ? nothing : push!(acclist, :initialized)
  for (clid, cliq) in treel.cliques
    if !(getCliqStatus(cliq) in acclist)
      return false
    end
  end
  return true
end

function isTreeSolvedUp(treel::BayesTree)::Bool
  for (clid, cliq) in treel.cliques
    if getCliqStatus(cliq) != :upsolved
      return false
    end
  end
  return true
end


"""
    $SIGNATURES

Return `::Bool` on whether all variables in this `cliq` are marginalzed.
"""
function isCliqMarginalizedFromVars(subfg::FactorGraph, cliq::Graphs.ExVertex)
  for vert in getCliqVars(subfg, cliq)
    if !isMarginalized(vert)
      return false
    end
  end
  return true
end

"""
    $SIGNATURES

Set the marginalized status of a clique.
"""
function setCliqAsMarginalized!(cliq::Graphs.ExVertex, status::Bool)
  if status
    getData(cliq).initialized = :marginalized
  else
    if getData(cliq).initialized == :marginalized
      @info "Reverting clique $(cliq.index) to assumed :downsolved status"
      getData(cliq).initialized = :downsolved
    else
      error("Unknown clique de-marginalization requist for clique $(cliq.index), current status: $(cliq.initialized)")
    end
  end
end

"""
    $SIGNATURES

Run through entire tree and set cliques as marginalized if all clique variables are marginalized.

Notes:
- TODO can be made fully parallel, consider converting for use with `@threads` `for`.
"""
function updateTreeCliquesAsMarginalizedFromVars!(fgl::FactorGraph, tree::BayesTree)::Nothing
  for (clid, cliq) in tree.cliques
    if isCliqMarginalizedFromVars(fgl, cliq)
      setCliqAsMarginalized!(cliq, true)
    end
  end
  nothing
end

"""
    $SIGNATURES

Reset the Bayes (Junction) tree so that a new upsolve can be performed.

Notes
- Will change previous clique status from `:downsolved` to `:initialized` only.
- Sets the color of tree clique to `lightgreen`.
"""
function resetTreeCliquesForUpSolve!(treel::BayesTree)::Nothing
  acclist = Symbol[:downsolved;]
  for (clid, cliq) in treel.cliques
    if getCliqStatus(cliq) in acclist
      setCliqStatus!(cliq, :initialized)
      setCliqDrawColor(cliq, "sienna")
    end
  end
  nothing
end

"""
    $SIGNATURES

Special internal function to try return the clique data if succesfully identified in `othertree::BayesTree`,
based on contents of `seeksSimilar::BayesTreeNodeData`.

Notes
- Used to identify and skip similar cliques (i.e. recycle computations)
"""
function attemptTreeSimilarClique(othertree::BayesTree, seeksSimilar::BayesTreeNodeData)::Graphs.ExVertex
  # inner convenience function for returning empty clique
  function EMPTYCLIQ()
    clq = ExVertex(-1,"null")
    clq.attributes["label"] = ""
    setData!(clq, emptyBTNodeData())
    return clq
  end

  # does the other clique even exist?
  seekFrontals = getCliqFrontalVarIds(seeksSimilar)
  if !hasCliq(othertree, seekFrontals[1])
    return EMPTYCLIQ()
  end

  # do the cliques share the same frontals?
  otherCliq = whichCliq(othertree, seekFrontals[1])
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



function tryCliqStateMachineSolve!(dfg::G,
                                   treel::BayesTree,
                                   i::Int;
                                   # cliqHistories;
                                   N::Int=100,
                                   oldtree::BayesTree=emptyBayesTree(),
                                   drawtree::Bool=false,
                                   limititers::Int=-1,
                                   downsolve::Bool=false,
                                   incremental::Bool=false,
                                   delaycliqs::Vector{Symbol}=Symbol[],
                                   recordcliqs::Vector{Symbol}=Symbol[]) where G <: AbstractDFG
  #
  clst = :na
  cliq = treel.cliques[i]
  syms = getCliqFrontalVarIds(cliq) # ids =
  oldcliq = attemptTreeSimilarClique(oldtree, getData(cliq))
  oldcliqdata = getData(oldcliq)
  opts = getSolverParams(dfg)
  # Base.rm(joinpath(opts.logpath,"logs/cliq$i"), recursive=true, force=true)
  mkpath(joinpath(opts.logpath,"logs/cliq$i/"))
  logger = SimpleLogger(open(joinpath(opts.logpath,"logs/cliq$i/log.txt"), "w+")) # NullLogger()
  # global_logger(logger)
  history = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}()
  recordthiscliq = length(intersect(recordcliqs,syms)) > 0
  delaythiscliq = length(intersect(delaycliqs,syms)) > 0
  try
    history = cliqInitSolveUpByStateMachine!(dfg, treel, cliq, N=N, drawtree=drawtree,
                                             oldcliqdata=oldcliqdata,
                                             limititers=limititers, downsolve=downsolve, recordhistory=recordthiscliq, incremental=incremental, delay=delaythiscliq, logger=logger )
    #
    # cliqHistories[i] = history
    if length(history) >= limititers && limititers != -1
      # @warn "writing logs/cliq$i/csm.txt"
      # @save "/tmp/cliqHistories/cliq$i.jld2" history
      fid = open(joinpath(opts.logpath,"logs/cliq$i/csm.txt"), "w")
      printCliqHistorySummary(fid, history)
      close(fid)
    end
    flush(logger.stream)
    close(logger.stream)
    # clst = getCliqStatus(cliq)
    # clst = cliqInitSolveUp!(dfg, treel, cliq, drawtree=drawtree, limititers=limititers )
  catch err
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

After solving, clique histories can be inserted back into the tree for later reference.
This function helps do the required assigment task.
"""
function assignTreeHistory!(treel::BayesTree, cliqHistories::Dict)
  for i in 1:length(treel.cliques)
    if haskey(cliqHistories, i)
      hist = cliqHistories[i]
      for i in 1:length(hist)
        hist[i][4].logger = SimpleLogger(stdout)
      end
      getData(treel.cliques[i]).statehistory=hist
    end
  end
end

function fetchCliqTaskHistoryAll!(smt, hist)
  for i in 1:length(smt)
    # sm = smt[i]
    hist[i] = fetch(smt[i])
  end
end

function fetchAssignTaskHistoryAll!(tree::BayesTree, smt)
  hist = Dict{Int, Vector{Tuple{DateTime,Int,Function,CliqStateMachineContainer}}}()
  fetchCliqTaskHistoryAll!(smt, hist)
  assignTreeHistory!(tree, hist)
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
function asyncTreeInferUp!(dfg::G,
                           treel::BayesTree;
                           oldtree::BayesTree=emptyBayesTree(),
                           drawtree::Bool=false,
                           N::Int=100,
                           limititers::Int=-1,
                           downsolve::Bool=false,
                           incremental::Bool=false,
                           skipcliqids::Vector{Symbol}=Symbol[],
                           delaycliqs::Vector{Symbol}=Symbol[],
                           recordcliqs::Vector{Symbol}=Symbol[] ) where G <: AbstractDFG
  #
  resetTreeCliquesForUpSolve!(treel)
  setTreeCliquesMarginalized!(dfg, treel)
  if drawtree
    pdfpath = joinpath(getSolverParams(dfg).logpath,"bt.pdf")
    drawTree(treel, show=false, filepath=pdfpath)
  end

  # queue all the tasks
  alltasks = Vector{Task}(undef, length(treel.cliques))
  # cliqHistories = Dict{Int,Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}()
  if !isTreeSolved(treel, skipinitialized=true)
    # @sync begin
      # duplicate int i into async (important for concurrency)
      for i in 1:length(treel.cliques)
        scsym = getCliqFrontalVarIds(treel.cliques[i])
        if length(intersect(scsym, skipcliqids)) == 0
          alltasks[i] = @async tryCliqStateMachineSolve!(dfg, treel, i, oldtree=oldtree, drawtree=drawtree, limititers=limititers, downsolve=downsolve, delaycliqs=delaycliqs, recordcliqs=recordcliqs, incremental=incremental, N=N)
        end # if
      end # for
    # end # sync
  end # if

  # post-hoc store possible state machine history in clique (without recursively saving earlier history inside state history)
  # assignTreeHistory!(treel, cliqHistories)

  # for i in 1:length(treel.cliques)
  #   if haskey(cliqHistories, i)
  #     getData(treel.cliques[i]).statehistory=cliqHistories[i]
  #   end
  # end

  return alltasks #, cliqHistories
end



"""
    $SIGNATURES

Perform tree based initialization of all variables not yet initialized in factor graph.

Related

asyncTreeInferUp!
"""
function initInferTreeUp!(dfg::G,
                          treel::BayesTree;
                          oldtree::BayesTree=emptyBayesTree(),
                          drawtree::Bool=false,
                          N::Int=100,
                          limititers::Int=-1,
                          downsolve::Bool=false,
                          incremental::Bool=false,
                          skipcliqids::Vector{Symbol}=Symbol[],
                          recordcliqs::Vector{Symbol}=Symbol[],
                          delaycliqs::Vector{Symbol}=Symbol[]) where G <: AbstractDFG
  #
  # revert :downsolved status to :initialized in preparation for new upsolve
  resetTreeCliquesForUpSolve!(treel)
  setTreeCliquesMarginalized!(dfg, treel)
  drawtree ? drawTree(treel, show=false, filepath=joinpath(getSolverParams(dfg).logpath,"bt.pdf")) : nothing

  # queue all the tasks
  alltasks = Vector{Task}(undef, length(treel.cliques))
  cliqHistories = Dict{Int,Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}}()
  if !isTreeSolved(treel, skipinitialized=true)
    @sync begin
      # duplicate int i into async (important for concurrency)
      for i in 1:length(treel.cliques)
        scsym = getCliqFrontalVarIds(treel.cliques[i])
        if length(intersect(scsym, skipcliqids)) == 0
          alltasks[i] = @async tryCliqStateMachineSolve!(dfg, treel, i, oldtree=oldtree, drawtree=drawtree, limititers=limititers, downsolve=downsolve, incremental=incremental, delaycliqs=delaycliqs, recordcliqs=recordcliqs) # N=N,
        end # if
      end # for
    end # sync
  end # if

  fetchCliqTaskHistoryAll!(alltasks, cliqHistories)

  # post-hoc store possible state machine history in clique (without recursively saving earlier history inside state history)
  assignTreeHistory!(treel, cliqHistories)
  # for i in 1:length(treel.cliques)
  #   if haskey(cliqHistories, i)
  #     hist = cliqHistories[i]
  #     for i in 1:length(hist)
  #       hist[i][4].logger = ConsoleLogger()
  #     end
  #     getData(treel.cliques[i]).statehistory=hist
  #   end
  # end

  return alltasks, cliqHistories
end
