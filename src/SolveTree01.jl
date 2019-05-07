


#global pidx
pidx = 1
pidl = 1
pidA = 1
thxl = nprocs() > 4 ? floor(Int,nprocs()*0.333333333) : 1

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
                                    wfac::Vector{<:AbstractString},
                                    vertid::Int,
                                    inmsgs::Array{NBPMessage,1},
                                    manis::T ) where {T <: Tuple}
  #
  for m in inmsgs
    for idx in keys(m.p)
      if idx == vertid
        pdi = m.p[vertid]
        push!(dens, manikde!(pdi.pts, pdi.bws, pdi.manifolds) ) # kde!(pdi.pts, pdi.bws)
        push!(wfac, "msg")
      end
      # TODO -- we can inprove speed of search for inner loop
    end
  end
  nothing
end

"""
    $(SIGNATURES)

Add all potentials associated with this clique and vertid to dens.
"""
function packFromLocalPotentials!(fgl::FactorGraph,
                                  dens::Vector{BallTreeDensity},
                                  wfac::Vector{AbstractString},
                                  cliq::Graphs.ExVertex,
                                  vertid::Int,
                                  N::Int,
                                  dbg::Bool=false )
  #
  for idfct in getData(cliq).potentials
    vert = getVert(fgl, idfct, api=localapi)
    data = getData(vert)
    # skip partials here, will be caught in packFromLocalPartials!
    if length( findall(data.fncargvID .== vertid) ) >= 1 && !data.fnc.partial
      p = findRelatedFromPotential(fgl, vert, vertid, N, dbg )
      push!(dens, p)
      push!(wfac, vert.label)
    end
  end
  nothing
end


function packFromLocalPartials!(fgl::FactorGraph,
      partials::Dict{Int, Vector{BallTreeDensity}},
      cliq::Graphs.ExVertex,
      vertid::Int,
      N::Int,
      dbg::Bool=false)
  #

  for idfct in getData(cliq).potentials
    vert = getVert(fgl, idfct, api=localapi)
    data = getData(vert)
    if length( findall(data.fncargvID .== vertid) ) >= 1 && data.fnc.partial
      p = findRelatedFromPotential(fgl, vert, vertid, N, dbg)
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
function productbelief(fg::FactorGraph,
                       vertid::Int,
                       dens::Vector{BallTreeDensity},
                       partials::Dict{Int, Vector{BallTreeDensity}},
                       N::Int;
                       dbg::Bool=false )
  #
  vert = getVert(fg, vertid, api=localapi)
  manis = getSofttype(vert).manifolds
  pGM = Array{Float64,2}(undef, 0,0)
  lennonp, lenpart = length(dens), length(partials)
  if lennonp > 1
    Ndims = Ndim(dens[1])
    @info "[$(lennonp)x$(lenpart)p,d$(Ndims),N$(N)],"
    pGM = prodmultiplefullpartials(dens, partials, Ndims, N, manis)
  elseif lennonp == 1 && lenpart >= 1
    Ndims = Ndim(dens[1])
    @info "[$(lennonp)x$(lenpart)p,d$(Ndims),N$(N)],"
    pGM = prodmultipleonefullpartials(dens, partials, Ndims, N, manis)
  elseif lennonp == 0 && lenpart >= 1
    denspts = getVal(fg,vertid,api=localapi)
    Ndims = size(denspts,1)
    @info "[$(lennonp)x$(lenpart)p,d$(Ndims),N$(N)],"
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
    @warn "Unknown density product on variable=$(vert.label), lennonp=$(lennonp), lenpart=$(lenpart)"
    pGM = Array{Float64,2}(undef, 0,1)
  end

  return pGM
end

function proposalbeliefs!(fgl::FactorGraph,
                          destvertid::Int,
                          factors::Vector{Graphs.ExVertex},
                          dens::Vector{BallTreeDensity},
                          partials::Dict{Int, Vector{BallTreeDensity}};
                          N::Int=100,
                          dbg::Bool=false,
                          api::DataLayerAPI=dlapi  )::Nothing
  #
  for fct in factors
    data = getData(fct)
    p = findRelatedFromPotential(fgl, fct, destvertid, N, dbg, api=api)
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
  end
  nothing
end

function predictbelief(fgl::FactorGraph,
                       destvert::ExVertex,
                       factors::Vector{Graphs.ExVertex};
                       N::Int=0,
                       dbg::Bool=false )
  #
  destvertid = destvert.index
  dens = Array{BallTreeDensity,1}()
  partials = Dict{Int, Vector{BallTreeDensity}}()

  # determine number of particles to draw from the marginal
  nn = N != 0 ? N : size(getVal(destvert),2)

  # get proposal beliefs
  proposalbeliefs!(fgl, destvertid, factors, dens, partials, N=nn, dbg=dbg)

  # take the product
  pGM = productbelief(fgl, destvertid, dens, partials, nn, dbg=dbg )

  return pGM
end

function predictbelief(fgl::FactorGraph,
                       destvertsym::Symbol,
                       factorsyms::Vector{Symbol};
                       N::Int=0,
                       dbg::Bool=false,
                       api::DataLayerAPI=IncrementalInference.localapi  )
  #
  factors = Graphs.ExVertex[]
  for fsym in factorsyms
    push!(factors, getVert(fgl, fgl.fIDs[fsym], api=api))
  end
  vert = getVert(fgl, destvertsym, api=api)

  # determine the number of particles to draw from the marginal
  nn = N != 0 ? N : size(getVal(vert),2)

  # do the belief prediction
  predictbelief(fgl, vert, factors, N=nn, dbg=dbg)
end

function predictbelief(fgl::FactorGraph,
                       destvertsym::Symbol,
                       factorsyms::Colon;
                       N::Int=0,
                       dbg::Bool=false,
                       api::DataLayerAPI=IncrementalInference.localapi  )
  #
  predictbelief(fgl, destvertsym, ls(fgl, destvertsym, api=api), N=N, api=api, dbg=dbg )
end

"""
    $(SIGNATURES)

Using factor graph object `fg`, project belief through connected factors
(convolution with conditional) to variable `sym` followed by a approximate functional product.

Return: product belief, full proposals, partial dimension proposals, labels
"""
function localProduct(fgl::FactorGraph,
                      sym::Symbol;
                      N::Int=100,
                      dbg::Bool=false,
                      api::DataLayerAPI=IncrementalInference.dlapi  )
  # TODO -- converge this function with predictbelief for this node
  # TODO -- update to use getVertId
  destvertid = fgl.IDs[sym] #destvert.index
  dens = Array{BallTreeDensity,1}()
  partials = Dict{Int, Vector{BallTreeDensity}}()
  lb = String[]
  fcts = Vector{Graphs.ExVertex}()
  cf = ls(fgl, sym, api=api)
  for f in cf
    vert = getVert(fgl, f, nt=:fnc, api=api)
    push!(fcts, vert)
    push!(lb, vert.label)
  end

  # get proposal beliefs
  proposalbeliefs!(fgl, destvertid, fcts, dens, partials, N=N, dbg=dbg, api=api)

  # take the product
  pGM = productbelief(fgl, destvertid, dens, partials, N, dbg=dbg )
  vert = getVert(fgl, sym, api=api)
  pp = AMP.manikde!(pGM, getSofttype(vert).manifolds )

  return pp, dens, partials, lb
end
localProduct(fgl::FactorGraph, lbl::T; N::Int=100, dbg::Bool=false) where {T <: AbstractString} = localProduct(fgl, Symbol(lbl), N=N, dbg=dbg)


"""
    $(SIGNATURES)

Initialize the belief of a variable node in the factor graph struct.
"""
function initVariable!(fgl::FactorGraph,
        sym::Symbol;
        N::Int=100,
        api::DataLayerAPI=IncrementalInference.dlapi )
  #

  vert = getVert(fgl, sym, api=api)
  # TODO -- this localapi is inconsistent, but get internal error due to problem with ls(fg, api=dlapi)
  belief,b,c,d  = localProduct(fgl, sym, api=localapi)
  pts = getPoints(belief)
  @show "initializing", sym, size(pts), Statistics.mean(pts,dims=2), Statistics.std(pts,dims=2)
  setVal!(vert, pts)
  api.updatevertex!(fgl, vert)

  nothing
end
function initializeNode!(fgl::FactorGraph,
                         sym::Symbol;
                         N::Int=100,
                         api::DataLayerAPI=IncrementalInference.dlapi )
  #
  @warn "initializeNode! has been deprecated in favor of initVariable!"
  initVariable!(fgl,sym,N=N,api=api )
end


"""
    $(SIGNATURES)

Perform one step of the minibatch clique Gibbs operation for solving the Chapman-Kolmogov
trasit integral -- here involving separate approximate functional convolution and
product operations.
"""
function cliqGibbs(fg::FactorGraph,
                   cliq::Graphs.ExVertex,
                   vertid::Int,
                   inmsgs::Array{NBPMessage,1},
                   N::Int,
                   dbg::Bool,
                   manis::T  ) where {T <: Tuple}
  #
  # several optimizations can be performed in this function TODO

  # consolidate NBPMessages and potentials
  dens = Array{BallTreeDensity,1}()
  partials = Dict{Int, Vector{BallTreeDensity}}()
  wfac = Vector{AbstractString}()
  packFromIncomingDensities!(dens, wfac, vertid, inmsgs, manis)
  packFromLocalPotentials!(fg, dens, wfac, cliq, vertid, N)
  packFromLocalPartials!(fg, partials, cliq, vertid, N, dbg)

  potprod = !dbg ? nothing : PotProd(vertid, getVal(fg,vertid,api=localapi), Array{Float64}(undef, 0,0), dens, wfac)
  pGM = productbelief(fg, vertid, dens, partials, N, dbg=dbg )
  if dbg  potprod.product = pGM  end

  # @info " "
  return pGM, potprod
end

"""
    $(SIGNATURES)

Iterate successive approximations of clique marginal beliefs by means
of the stipulated proposal convolutions and products of the functional objects
for tree clique `cliq`.
"""
function fmcmc!(fgl::FactorGraph,
                cliq::Graphs.ExVertex,
                fmsgs::Vector{NBPMessage},
                IDs::Vector{Int},
                N::Int,
                MCMCIter::Int,
                dbg::Bool=false,
                api::DataLayerAPI=dlapi )
  #
    @info "---------- successive fnc approx ------------$(cliq.attributes["label"])"
    # repeat several iterations of functional Gibbs sampling for fixed point convergence
    if length(IDs) == 1
        MCMCIter=1
    end
    mcmcdbg = Array{CliqGibbsMC,1}()

    for iter in 1:MCMCIter
      # iterate through each of the variables, KL-divergence tolerence would be nice test here
      @info "#$(iter)\t -- "
      dbgvals = !dbg ? nothing : CliqGibbsMC([], Symbol[])
      # @show IDs
      for vertid in IDs
        vert = getVert(fgl, vertid, api=api)
        # @show vert.index, vert.label
        if !getData(vert).ismargin
          # we'd like to do this more pre-emptive and then just execute -- just point and skip up only msgs
          densPts, potprod = cliqGibbs(fgl, cliq, vertid, fmsgs, N, dbg, getSofttype(vert).manifolds) #cliqGibbs(fg, cliq, vertid, fmsgs, N)
          if size(densPts,1)>0
            updvert = getVert(fgl, vertid, api=dlapi)  # TODO --  can we remove this duplicate getVert?
            setValKDE!(updvert, densPts)
            # Go update the datalayer TODO -- excessive for general case, could use local and update remote at end
            dlapi.updatevertex!(fgl, updvert)
            # fgl.v[vertid].attributes["val"] = densPts
            if dbg
              push!(dbgvals.prods, potprod)
              push!(dbgvals.lbls, Symbol(updvert.label))
            end
          end
        end
      end
      !dbg ? nothing : push!(mcmcdbg, dbgvals)
      # @info ""
    end

    # populate dictionary for return NBPMessage in multiple dispatch
    # TODO -- change to EasyMessage dict
    d = Dict{Int,EasyMessage}() # Array{Float64,2}
    for vertid in IDs
      # TODO reduce to local fg only
      vert = getVert(fgl,vertid, api=api)
      pden = getKDE(vert)
      bws = vec(getBW(pden)[:,1])
      manis = getSofttype(vert).manifolds
      d[vertid] = EasyMessage(getVal(vert), bws, manis)
      # d[vertid] = getVal(dlapi.getvertex(fgl,vertid)) # fgl.v[vertid]
    end
    @info "fmcmc! -- finished on $(cliq.attributes["label"])"

    return mcmcdbg, d
end

function upPrepOutMsg!(d::Dict{Int,EasyMessage}, IDs::Array{Int,1}) #Array{Float64,2}
  @info "Outgoing msg density on: "
  len = length(IDs)
  m = NBPMessage(Dict{Int,EasyMessage}())
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
function treeProductUp(fg::FactorGraph,
                       tree::BayesTree,
                       cliq::Symbol,
                       sym::Symbol;
                       N::Int=100,
                       dbg::Bool=false  )
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
  pGM, potprod = cliqGibbs( fg, cliq, vertid, upmsgssym, N, dbg, manis )

  return pGM, potprod
end


"""
    $(SIGNATURES)

Calculate a fresh---single step---approximation to the variable `sym` in clique `cliq` as though during the downward message passing.  The full inference algorithm may repeatedly calculate successive apprimxations to the variable based on the structure of variables, factors, and incoming messages to this clique.
Which clique to be used is defined by frontal variable symbols (`cliq` in this case) -- see `whichCliq(...)` for more details.  The `sym` symbol indicates which symbol of this clique to be calculated.  **Note** that the `sym` variable must appear in the clique where `cliq` is a frontal variable.
"""
function treeProductDwn(fg::FactorGraph,
                        tree::BayesTree,
                        cliq::Symbol,
                        sym::Symbol;
                        N::Int=100,
                        dbg::Bool=false  )
  #
  cliq = whichCliq(tree, cliq)
  cliqdata = getData(cliq)

  # get the local variable id::Int identifier
  vertid = fg.IDs[sym]

  # get all the incoming (upward) messages from the tree cliques
  # convert incoming messages to Int indexed format (semi-legacy format)
  cl = parentCliq(tree, cliq)
  msgdict = dwnMsg(cl[1])
  dict = Dict{Int, EasyMessage}()
  for (dsy, btd) in msgdict
      dict[fg.IDs[dsy]] = convert(EasyMessage, btd)
  end
  dwnmsgssym = NBPMessage[NBPMessage(dict);]

  # perform the actual computation
  pGM, potprod = cliqGibbs( fg, cliq, vertid, dwnmsgssym, N, dbg )

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
function upGibbsCliqueDensity(inp::ExploreTreeType{T},
                              N::Int=100,
                              dbg::Bool=false,
                              iters::Int=3) where {T}
  #
  @info "up w $(length(inp.sendmsgs)) msgs"
  # Local mcmc over belief functions
  # this is so slow! TODO Can be ignored once we have partial working
  # loclfg = nprocs() < 2 ? deepcopy(inp.fg) : inp.fg

  # TODO -- some weirdness with: d,. = d = ., nothing
  mcmcdbg = Array{CliqGibbsMC,1}()
  d = Dict{Int,EasyMessage}()

  priorprods = Vector{CliqGibbsMC}()

  cliqdata = getData(inp.cliq)

  # use nested structure for more fficient Chapman-Kolmogorov solution approximation
  if false
    IDS = [cliqdata.frontalIDs;cliqdata.conditIDs] #inp.cliq.attributes["frontalIDs"]
    mcmcdbg, d = fmcmc!(inp.fg, inp.cliq, inp.sendmsgs, IDS, N, iters, dbg)
  else
    dummy, d = fmcmc!(inp.fg, inp.cliq, inp.sendmsgs, cliqdata.directFrtlMsgIDs, N, 1)
    if length(cliqdata.msgskipIDs) > 0
      dummy, dd = fmcmc!(inp.fg, inp.cliq, inp.sendmsgs, cliqdata.msgskipIDs, N, 1)
      for md in dd d[md[1]] = md[2]; end
    end
    # NOTE -- previous mistake, must iterate over directsvarIDs also
    if length(cliqdata.itervarIDs) > 0
      mcmcdbg, ddd = fmcmc!(inp.fg, inp.cliq, inp.sendmsgs, cliqdata.itervarIDs, N, iters, dbg)
      for md in ddd d[md[1]] = md[2]; end
    end
    if length(cliqdata.directPriorMsgIDs) > 0
      doids = setdiff(cliqdata.directPriorMsgIDs, cliqdata.msgskipIDs)
      priorprods, dddd = fmcmc!(inp.fg, inp.cliq, inp.sendmsgs, doids, N, 1, dbg)
      for md in dddd d[md[1]] = md[2]; end
    end
  end

  #m = upPrepOutMsg!(inp.fg, inp.cliq, inp.sendmsgs, condids, N)
  m = upPrepOutMsg!(d, cliqdata.conditIDs)

  outmsglbl = Dict{Symbol, Int}()
  if dbg
    for (ke, va) in m.p
      outmsglbl[Symbol(inp.fg.g.vertices[ke].label)] = ke
    end
  end

  # prepare and convert upward belief messages
  upmsgs = Dict{Symbol, BallTreeDensity}()
  for (ke, va) in m.p
    msgsym = Symbol(inp.fg.g.vertices[ke].label)
    upmsgs[msgsym] = convert(BallTreeDensity, va)
  end
  setUpMsg!(inp.cliq, upmsgs)

  # flag cliq as definitely being initialized
  cliqdata.upsolved = true

  mdbg = !dbg ? DebugCliqMCMC() : DebugCliqMCMC(mcmcdbg, m, outmsglbl, priorprods)
  return UpReturnBPType(m, mdbg, d, upmsgs, true)
end



function dwnPrepOutMsg(fg::FactorGraph, cliq::Graphs.ExVertex, dwnMsgs::Array{NBPMessage,1}, d::Dict{Int, EasyMessage}) #Array{Float64,2}
    # pack all downcoming conditionals in a dictionary too.
    if cliq.index != 1
      @info "Dwn msg keys $(keys(dwnMsgs[1].p))"
    end # ignore root, now incoming dwn msg
    @info "Outgoing msg density on: "
    m = NBPMessage(Dict{Int,EasyMessage}())
    i = 0
    for vid in cliq.attributes["data"].frontalIDs
      m.p[vid] = deepcopy(d[vid]) # TODO -- not sure if deepcopy is required
    end
    for cvid in cliq.attributes["data"].conditIDs
        i+=1
        # TODO -- convert to points only since kde replace by rkhs in future
        m.p[cvid] = deepcopy(dwnMsgs[1].p[cvid]) # TODO -- maybe this can just be a union(,)
    end
    return m
end

"""
    $SIGNATURES

Perform Chapman-Kolmogorov transit integral approximation for `cliq` in downward pass direction.

Note
- Only update frontal variables of the clique.
"""
function downGibbsCliqueDensity(fg::FactorGraph,
                                cliq::Graphs.ExVertex,
                                dwnMsgs::Array{NBPMessage,1},
                                N::Int=100,
                                MCMCIter::Int=3,
                                dbg::Bool=false  )
    #
    # TODO standardize function call to have similar stride to upGibbsCliqueDensity
    @info "down"
    mcmcdbg, d = fmcmc!(fg, cliq, dwnMsgs, getFrontals(cliq), N, MCMCIter, dbg)
    m = dwnPrepOutMsg(fg, cliq, dwnMsgs, d)

    outmsglbl = Dict{Symbol, Int}()
    if dbg
      for (ke, va) in m.p
        outmsglbl[Symbol(fg.g.vertices[ke].label)] = ke
      end
    end

    # Always keep dwn messages in cliq data
    dwnkeepmsgs = Dict{Symbol, BallTreeDensity}()
    for (ke, va) in m.p
      msgsym = Symbol(fg.g.vertices[ke].label)
      dwnkeepmsgs[msgsym] = convert(BallTreeDensity, va)
    end
    setDwnMsg!(cliq, dwnkeepmsgs)

    # down solving complete, set flag
    getData(cliq).downsolved = true

    mdbg = !dbg ? DebugCliqMCMC() : DebugCliqMCMC(mcmcdbg, m, outmsglbl, CliqGibbsMC[])
    return DownReturnBPType(m, mdbg, d, dwnkeepmsgs)
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
function updateFGBT!(fg::FactorGraph,
                     bt::BayesTree,
                     cliqID::Int,
                     ddt::DownReturnBPType;
                     dbg::Bool=false,
                     fillcolor::String=""  )
    #
    cliq = bt.cliques[cliqID]
    if dbg
      cliq.attributes["debugDwn"] = deepcopy(ddt.dbgDwn)
    end
    setDwnMsg!(cliq, ddt.keepdwnmsgs)
    # TODO move to drawTree
    if fillcolor != ""
      setCliqDrawColor(cliq, fillcolor)
    end
    for dat in ddt.IDvals
      #TODO -- should become an update call
        updvert = dlapi.getvertex(fg,dat[1])
        setValKDE!(updvert, deepcopy(dat[2])) # TODO -- not sure if deepcopy is required
        dlapi.updatevertex!(fg, updvert, updateMAPest=true)
    end
    nothing
end

"""
    $(SIGNATURES)

Update cliq `cliqID` in Bayes (Juction) tree `bt` according to contents of `urt` -- intended use is to update main clique after a upward belief propagation computation has been completed per clique.
"""
function updateFGBT!(fg::FactorGraph,
                     bt::BayesTree,
                     cliqID::Int,
                     urt::UpReturnBPType;
                     dbg::Bool=false, fillcolor::String="" )
  #
  cliq = bt.cliques[cliqID]
  cliq = bt.cliques[cliqID]
  if dbg
    cliq.attributes["debug"] = deepcopy(urt.dbgUp)
  end
  setUpMsg!(cliq, urt.keepupmsgs)
  # move to drawTree
  if fillcolor != ""
    setCliqDrawColor(cliq, fillcolor)
  end
  for dat in urt.IDvals
    updvert = dlapi.getvertex(fg,dat[1])
    setValKDE!(updvert, deepcopy(dat[2])) # (fg.v[dat[1]], ## TODO -- not sure if deepcopy is required
    dlapi.updatevertex!(fg, updvert, updateMAPest=true)
  end
  @info "updateFGBT! up -- finished updating $(cliq.attributes["label"])"
  nothing
end

"""
    $(SIGNATURES)

Pass NBPMessages back down the tree -- pre order tree traversal.
"""
function downMsgPassingRecursive(inp::ExploreTreeType{T}; N::Int=100, dbg::Bool=false, drawpdf::Bool=false) where {T}
  @info "====================== Clique $(inp.cliq.attributes["label"]) ============================="

  mcmciter = inp.prnt != nothing ? 3 : 0; # skip mcmc in root on dwn pass
  rDDT = downGibbsCliqueDensity(inp.fg, inp.cliq, inp.sendmsgs, N, mcmciter, dbg) #dwnMsg
  updateFGBT!(inp.fg, inp.bt, inp.cliq.index, rDDT, dbg=dbg, fillcolor="lightblue")
  setCliqStatus!(inp.cliq, :downsolved)
  drawpdf ? drawTree(inp.bt) : nothing

  # rr = Array{Future,1}()
  # pcs = procs()

  ddt=nothing
  for child in out_neighbors(inp.cliq, inp.bt.bt)
    ett = ExploreTreeType(inp.fg, inp.bt, child, inp.cliq, [rDDT.dwnMsg])#inp.fg
    ddt = downMsgPassingRecursive( ett , N=N, dbg=dbg )
    drawpdf ? drawTree(inp.bt) : nothing
  end

  # return modifications to factorgraph to calling process
  return ddt
end

"""
    $SIGNATURES

Get and return upward belief messages as stored in child cliques from `treel::BayesTree`.
"""
function getCliqChildMsgsUp(fg_::FactorGraph, treel::BayesTree, cliq::Graphs.ExVertex, ::Type{EasyMessage})
  childmsgs = NBPMessage[]
  for child in getChildren(treel, cliq)
    nbpchild = NBPMessage(Dict{Int,EasyMessage}())
    for (key, bel) in getUpMsgs(child)
      @info "$(current_task()) Clique $(cliq.index), child cliq $(child.index), getCliqChildMsgsUp -- key=$(key)"
      id = fg_.IDs[key]
      manis = getManifolds(fg_, id)
      nbpchild.p[id] = convert(EasyMessage, bel, manis)
    end
    push!(childmsgs, nbpchild)
  end
  return childmsgs
end


"""
    $SIGNATURES

Approximate Chapman-Kolmogorov transit integral and return separator marginals as messages to pass up the Bayes (Junction) tree, along with additional clique operation values for debugging.

Notes
-----
- `onduplicate=true` by default internally uses deepcopy of factor graph and Bayes tree, and does **not** update the given objects.  Set false to update `fgl` and `treel` during compute.
"""
function approxCliqMarginalUp!(fgl::FactorGraph,
                               treel::BayesTree,
                               csym::Symbol,
                               onduplicate=true;
                               N::Int=100,
                               dbg::Bool=false,
                               iters::Int=3,
                               drawpdf::Bool=false  )
  #
  fg_ = onduplicate ? deepcopy(fgl) : fgl
  onduplicate ? (@warn "rebuilding new Bayes tree on deepcopy of factor graph") : nothing
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

  @info "=== start Clique $(cliq.attributes["label"]) ======================"
  ett = ExploreTreeType(fg_, tree_, cliq, nothing, childmsgs)
  urt = upGibbsCliqueDensity(ett, N, dbg, iters)
  updateFGBT!(ett.fg, ett.bt, ett.cliq.index, urt, dbg=dbg, fillcolor="pink")
  drawpdf ? drawTree(tree_) : nothing
  @info "=== end Clique $(cliq.attributes["label"]) ========================"
  urt
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
                            drawpdf::Bool=false   )
  #
  approxCliqMarginalUp!(fgl, treel, csym, onduplicate; N=N, dbg=dbg, iters=iters, drawpdf=drawpdf )
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


# post order tree traversal and build potential functions
function upMsgPassingRecursive(inp::ExploreTreeType{T}; N::Int=100, dbg::Bool=false, drawpdf::Bool=false) where {T}
    @info "Start Clique $(inp.cliq.attributes["label"]) ============================="
    childMsgs = Array{NBPMessage,1}()

    outnei = out_neighbors(inp.cliq, inp.bt.bt)
    len = length(outnei)
    for child in outnei
        ett = ExploreTreeType(inp.fg, inp.bt, child, inp.cliq, NBPMessage[])
        @info "upMsgRec -- calling new recursive on $(ett.cliq.attributes["label"])"
        newmsgs = upMsgPassingRecursive(  ett, N=N, dbg=dbg ) # newmsgs
        @info "upMsgRec -- finished with $(ett.cliq.attributes["label"]), w $(keys(newmsgs.p)))"
        push!(  childMsgs, newmsgs )
    end

    @info "====================== Clique $(inp.cliq.attributes["label"]) ============================="
    ett = ExploreTreeType(inp.fg, inp.bt, inp.cliq, nothing, childMsgs)

    urt = upGibbsCliqueDensity(ett, N, dbg) # upmsgdict
    updateFGBT!(inp.fg, inp.bt, inp.cliq.index, urt, dbg=dbg, fillcolor="lightblue")
    drawpdf ? drawTree(inp.bt) : nothing
    @info "End Clique $(inp.cliq.attributes["label"]) ============================="
    urt.upMsgs
end


function downGibbsCliqueDensity(inp::ExploreTreeType{T},
                                N::Int=100,
                                dbg::Bool=false  ) where {T}
  #
  @info "=================== Iter Clique $(inp.cliq.attributes["label"]) ==========================="
  mcmciter = inp.prnt != nothing ? 3 : 0
  return downGibbsCliqueDensity(inp.fg, inp.cliq, inp.sendmsgs, N, mcmciter, dbg)
end

function prepDwnPreOrderStack!(bt::BayesTree, parentStack::Array{Graphs.ExVertex,1})
  # dwn message passing function
  nodedata = nothing
  tempStack = Array{Graphs.ExVertex,1}()
  push!(tempStack, parentStack[1])

  while ( length(tempStack) != 0 || nodedata != nothing )
    if nodedata != nothing
      for child in out_neighbors(nodedata, bt.bt) #nodedata.cliq, nodedata.bt.bt
          ## TODO -- don't yet have the NBPMessage results to propagate belief
          ## should only construct ett during processing
          push!(parentStack, child) #ett
          push!(tempStack, child) #ett
      end
      nodedata = nothing # its over but we still need to complete computation in each of the leaves of the tree
    else
      nodedata = tempStack[1]
      deleteat!(tempStack, 1)
    end
  end
  nothing
end

function findVertsAssocCliq(fgl::FactorGraph, cliq::Graphs.ExVertex)

  cliqdata = getData(cliq)
  IDS = [cliqdata.frontalIDs; cliqdata.conditIDs] #inp.cliq.attributes["frontalIDs"]


  @error "findVertsAssocCliq -- not completed yet"
  nothing
end

function partialExploreTreeType(pfg::FactorGraph, pbt::BayesTree, cliqCursor::Graphs.ExVertex, prnt, pmsgs::Array{NBPMessage,1})
    # info("starting pett")
    # TODO -- expand this to grab only partial subsection from the fg and bt data structures


    if length(pmsgs) < 1
      return ExploreTreeType(pfg, pbt, cliqCursor, prnt, NBPMessage[])
    else
      return ExploreTreeType(pfg, pbt, cliqCursor, prnt, pmsgs)
    end
    nothing
end

function dispatchNewDwnProc!(fg::FactorGraph,
                             bt::BayesTree,
                             parentStack::Array{Graphs.ExVertex,1},
                             stkcnt::Int,
                             refdict::Dict{Int,Future};
                             N::Int=100,
                             dbg::Bool=false,
                             drawpdf::Bool=false  )
  #
  cliq = parentStack[stkcnt]
  while !haskey(refdict, cliq.index) # nodedata.cliq
    sleep(0.25)
  end

  rDDT = fetch(refdict[cliq.index]) #nodedata.cliq
  delete!(refdict, cliq.index) # nodedata

  if rDDT != nothing
    updateFGBT!(fg, bt, cliq.index, rDDT, dbg=dbg, fillcolor="lightblue")
    setCliqStatus!(cliq, :downsolved)
    drawpdf ? drawTree(bt) : nothing
  end

  emptr = BayesTree(nothing, 0, Dict{Int,Graphs.ExVertex}(), Dict{String,Int}());

  for child in out_neighbors(cliq, bt.bt) # nodedata.cliq, nodedata.bt.bt
      haskey(refdict, child.index) ? error("dispatchNewDwnProc! -- why you already have dwnremoteref?") : nothing
      ett = partialExploreTreeType(fg, emptr, child, cliq, [rDDT.dwnMsg]) # bt
      refdict[child.index] = remotecall(downGibbsCliqueDensity, upp2() , ett, N) # Julia 0.5 swapped order
  end
  nothing
end

"""
    $SIGNATURES

Downward message passing on Bayes (Junction) tree.

Notes
- Simultaenously launches as many async dispatches to remote processes as there are cliques in the tree.
"""
function processPreOrderStack!(fg::FactorGraph,
                               bt::BayesTree,
                               parentStack::Array{Graphs.ExVertex,1},
                               refdict::Dict{Int,Future};
                               N::Int=100,
                               dbg::Bool=false,
                               drawpdf::Bool=false )
  #
    # dwn message passing function for iterative tree exploration
    stkcnt = 0

    @sync begin
      sendcnt = 1:length(parentStack) # separate memory for remote calls
      for i in 1:sendcnt[end]
          @async try
            dispatchNewDwnProc!(fg, bt, parentStack, sendcnt[i], refdict, N=N, dbg=dbg, drawpdf=drawpdf) # stkcnt ##pidxI,nodedata
          catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)
          end
      end
    end
    nothing
end

function downMsgPassingIterative!(startett::ExploreTreeType{T};
                                  N::Int=100,
                                  dbg::Bool=false,
                                  drawpdf::Bool=false  ) where {T}
  #
  # this is where we launch the downward iteration process from
  parentStack = Array{Graphs.ExVertex,1}()
  refdict = Dict{Int,Future}()

  # start at the given clique in the tree -- shouldn't have to be the root.
  pett = partialExploreTreeType(startett.fg, startett.bt, startett.cliq,
                                        startett.prnt, startett.sendmsgs)
  refdict[startett.cliq.index] = remotecall(downGibbsCliqueDensity, upp2(), pett, N)  # for Julia 0.5

  push!(parentStack, startett.cliq )

  prepDwnPreOrderStack!(startett.bt, parentStack)
  processPreOrderStack!(startett.fg, startett.bt, parentStack, refdict, N=N, dbg=dbg, drawpdf=drawpdf )

  @info "dwnward leftovers, $(keys(refdict))"
  nothing
end

function prepPostOrderUpPassStacks!(bt::BayesTree,
                                    parentStack::Array{Graphs.ExVertex,1},
                                    childStack::Array{Graphs.ExVertex,1}  )
  # upward message passing preparation
  while ( length(parentStack) != 0 )
      #2.1 Pop a node from first stack and push it to second stack
      cliq = parentStack[end]
      deleteat!(parentStack, length(parentStack))
      push!(childStack, cliq )

      #2.2 Push left and right children of the popped node to first stack
      for child in out_neighbors(cliq, bt.bt) # nodedata.cliq, nodedata.bt.bt
          push!(parentStack, child )
      end
  end
  for child in childStack
    @show child.attributes["label"]
  end
  nothing
end

"""
    $SIGNATURES

Asynchronously perform up message passing, based on previoulsy prepared `chldstk::Vector{ExVertex}`.
"""
function asyncProcessPostStacks!(fgl::FactorGraph,
                                 bt::BayesTree,
                                 chldstk::Vector{Graphs.ExVertex},
                                 stkcnt::Int,
                                 refdict::Dict{Int,Future};
                                 N::Int=100,
                                 dbg::Bool=false,
                                 drawpdf::Bool=false  )
  #
  if stkcnt == 0
    @info "asyncProcessPostStacks! ERROR stkcnt=0"
    error("asyncProcessPostStacks! stkcnt=0")
  end
  cliq = chldstk[stkcnt]
  gomulti = true
  @info "Start Clique $(cliq.attributes["label"]) ============================="
  childMsgs = Array{NBPMessage,1}()
  ur = nothing
  for child in out_neighbors(cliq, bt.bt)
      @info "asyncProcessPostStacks -- $(stkcnt), cliq=$(cliq.attributes["label"]), start on child $(child.attributes["label"]) haskey=$(haskey(child.attributes, "remoteref"))"
        while !haskey(refdict, child.index)
          # info("Sleeping $(cliq.attributes["label"]) on lack of remoteref from $(child.attributes["label"])")
          # @show child.index, keys(refdict)
          sleep(0.25)
        end

      if gomulti
        ur = fetch(refdict[child.index])
      else
        ur = child.attributes["remoteref"]
      end
      updateFGBT!( fgl, bt, child.index, ur, dbg=dbg, fillcolor="pink" ) # deep copies happen in the update function
      drawpdf ? drawTree(bt) : nothing
      #delete!(child.attributes, "remoteref")

      push!(childMsgs, ur.upMsgs)

  end
  @info "====================== Clique $(cliq.attributes["label"]) ============================="
  emptr = BayesTree(nothing, 0, Dict{Int,Graphs.ExVertex}(), Dict{String,Int}());
  pett = partialExploreTreeType(fgl, emptr, cliq, nothing, childMsgs) # bt   # parent cliq pointer is not needed here, fix Graphs.jl first

  if haskey(cliq.attributes, "remoteref")
      @info "asyncProcessPostStacks! -- WHY YOU ALREADY HAVE REMOTEREF?"
  end

  newprocid = upp2()
  if gomulti
    refdict[cliq.index] = remotecall(upGibbsCliqueDensity, newprocid, pett, N, dbg ) # swap order for Julia 0.5
  else
    # bad way to do non multi test
    cliq.attributes["remoteref"] = upGibbsCliqueDensity(pett, N, dbg)
  end

  # delete as late as possible, but could happen sooner
  for child in out_neighbors(cliq, bt.bt)
      # delete!(child.attributes, "remoteref")
      delete!(refdict, child.index)
  end

  @info "End Clique $(cliq.attributes["label"]) ============================="
  nothing
end

"""
    $SIGNATURES

Multiprocess upward belief propagation message passing function, using async tasks.

Notes
- asyncs used to wrap remotecall for multicore.
- separate multithreaded calls can occur on each separate process.
"""
function processPostOrderStacks!(fg::FactorGraph,
                                 bt::BayesTree,
                                 childStack::Array{Graphs.ExVertex,1};
                                 N::Int=100,
                                 dbg::Bool=false,
                                 drawpdf::Bool=false  )
  #
  refdict = Dict{Int,Future}()

  stkcnt = length(childStack)
  @sync begin
    sendcnt = stkcnt:-1:1 # separate stable memory
    for i in 1:stkcnt
        @async asyncProcessPostStacks!(fg, bt, childStack, sendcnt[i], refdict, N=N, dbg=dbg, drawpdf=drawpdf ) # deepcopy(stkcnt)
    end
  end
  @info "processPostOrderStacks! -- THIS ONLY HAPPENS AFTER SYNC"
  # we still need to fetch the root node computational output
  if true
    ur = fetch(refdict[childStack[1].index])
  else
    ur = childStack[1].attributes["remoteref"]
  end
  # delete!(childStack[1].attributes, "remoteref") # childStack[1].cliq
  delete!(refdict, childStack[1].index)

  @info "upward leftovers, $(keys(refdict))"

  updateFGBT!(fg, bt, childStack[1].index, ur, dbg=dbg, fillcolor="pink" ) # nodedata
  drawpdf ? drawTree(bt) : nothing
  nothing
end


"""
    $SIGNATURES

Return clique pointers for the given order in which they will be solved (sequentially).
"""
function getCliqOrderUpSolve(treel::BayesTree, startcliq=treel.cliques[1])
  # http://www.geeksforgeeks.org/iterative-postorder-traversal/
  # this is where we launch the downward iteration process from
  parentStack = Vector{Graphs.ExVertex}()
  childStack = Vector{Graphs.ExVertex}()
  #Loop while first stack is not empty
  push!(parentStack, startcliq)
  # Starting at the root means we have a top down view of the tree
  prepPostOrderUpPassStacks!(treel, parentStack, childStack)
  return childStack
end

"""
    $SIGNATURES

Return clique pointers for the given order in which they will be solved (sequentially).
"""
getTreeCliqSolveOrderUp(treel::BayesTree, startcliq=treel.cliques[1]) = getCliqOrderUpSolve(treel, startcliq)

"""
    $SIGNATURES

Perform upward message passing (multi-process) algorithm for sum-product solution from leaves to root of the tree.

Notes:
* inspired by http://www.geeksforgeeks.org/iterative-postorder-traversal/
* this is where downward iteration process is launched from.
"""
function upMsgPassingIterative!(startett::ExploreTreeType{T};
                                N::Int=100,
                                dbg::Bool=false,
                                drawpdf::Bool=false  ) where {T}
  #
  childStack = getCliqOrderUpSolve(startett.bt, startett.cliq)
  # Starting at the root means we have a top down view of the tree
  processPostOrderStacks!(startett.fg, startett.bt, childStack, N=N, dbg=dbg, drawpdf=drawpdf)
  nothing
end
# for (ids, cliq) in treel.cliques
#   getData(cliq).initialized = :initialized
# end

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

Perform tree based initialization of all variables not yet initialized in factor graph.
"""
function initInferTreeUp!(fgl::FactorGraph,
                          treel::BayesTree;
                          drawtree::Bool=false,
                          N::Int=100,
                          limititers::Int=-1  )
  # revert :downsolved status to :initialized in preparation for new upsolve
  resetTreeCliquesForUpSolve!(treel)
  setTreeCliquesMarginalized!(fgl, treel)
  drawtree ? drawTree(treel, show=false) : nothing

  # queue all the tasks
  alltasks = Vector{Task}(undef, length(treel.cliques))
  @sync begin
    if !isTreeSolved(treel, skipinitialized=true)
      for i in 1:length(treel.cliques)
        alltasks[i] = @async begin
          clst = :na
          try
            cliq = treel.cliques[i]
            history = cliqInitSolveUpByStateMachine!(fgl, treel, cliq, drawtree=drawtree, limititers=limititers )
            clst = getCliqStatus(cliq)
            # clst = cliqInitSolveUp!(fgl, treel, cliq, drawtree=drawtree, limititers=limititers )
          catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)
            error(err)
          end
          if !(clst in [:upsolved; :downsolved; :marginalized])
            error("Clique $(cliq.index), initInferTreeUp! -- cliqInitSolveUp! did not arrive at the desired solution statu: $clst")
          end
        end # async
      end # for
    end # if
  end # sync

  return alltasks
end


"""
    $SIGNATURES

Perform up and down message passing (multi-process) algorithm for full sum-product solution of all continuous marginal beliefs.
"""
function inferOverTree!(fgl::FactorGraph,
                        bt::BayesTree;
                        N::Int=100,
                        dbg::Bool=false,
                        drawpdf::Bool=false,
                        treeinit::Bool=false  )::Nothing
  #
  @info "Batch rather than incremental solving over the Bayes (Junction) tree."
  setAllSolveFlags!(bt, false)
  @info "Ensure all nodes are initialized"
  inittasks = Vector{Task}()
  if treeinit
    @info "Do tree based init-inference on tree"
    inittasks = initInferTreeUp!(fgl, bt, N=N, drawtree=drawpdf)
    @info "Finished tree based upward init-inference"
  else
    ensureAllInitialized!(fgl)
  end
  if !treeinit # !isTreeSolvedUp(bt)
    @info "Do multi-process upward pass of inference on tree"
    upMsgPassingIterative!(ExploreTreeType(fgl, bt, bt.cliques[1], nothing, NBPMessage[]),N=N, dbg=dbg, drawpdf=drawpdf);
  end
  @info "Do multi-process downward pass of inference on tree"
  downMsgPassingIterative!(ExploreTreeType(fgl, bt, bt.cliques[1], nothing, NBPMessage[]),N=N, dbg=dbg, drawpdf=drawpdf);
  nothing
end

"""
    $SIGNATURES

Perform up and down message passing (single process, recursive) algorithm for full sum-product solution of all continuous marginal beliefs.
"""
function inferOverTreeR!(fgl::FactorGraph,
                         bt::BayesTree;
                         N::Int=100,
                         dbg::Bool=false,
                         drawpdf::Bool=false,
                         treeinit::Bool=false  )::Nothing
  #
  @info "Batch rather than incremental solving over the Bayes (Junction) tree."
  setAllSolveFlags!(bt, false)
  @info "Ensure all nodes are initialized"
  if treeinit
    inittasks = initInferTreeUp!(fgl, bt, N=N, drawtree=drawpdf)
  else
    @info "Do conventional recursive up inference over tree"
    ensureAllInitialized!(fgl)
    upMsgPassingRecursive(ExploreTreeType(fgl, bt, bt.cliques[1], nothing, NBPMessage[]), N=N, dbg=dbg, drawpdf=drawpdf);
  end
  @info "Do recursive down inference over tree"
  downMsgPassingRecursive(ExploreTreeType(fgl, bt, bt.cliques[1], nothing, NBPMessage[]), N=N, dbg=dbg, drawpdf=drawpdf);
  nothing
end
