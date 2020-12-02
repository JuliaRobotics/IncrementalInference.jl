


"""
    $(SIGNATURES)

Multiply different dimensions from partial constraints individually.
"""
function productpartials!(pGM::Array{Float64,2},
                          dummy::BallTreeDensity,
                          partials::Dict{Int, Vector{BallTreeDensity}},
                          manis::T  ) where {T <: Tuple}
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
function prodmultiplefullpartials(dens::Vector{BallTreeDensity},
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
- Incorporate ApproxManifoldProducts to process variables in individual batches.
"""
function productbelief( dfg::AbstractDFG,
                        vertlabel::Symbol,
                        dens::Vector{BallTreeDensity},
                        partials::Dict{Int, Vector{BallTreeDensity}},
                        N::Int;
                        dbg::Bool=false,
                        logger=ConsoleLogger()  )
  #
  vert = DFG.getVariable(dfg, vertlabel)
  manis = getVariableType(vert) |> getManifolds
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
function proposalbeliefs!(dfg::AbstractDFG,
                          destvertlabel::Symbol,
                          factors::AbstractVector{<:DFGFactor},
                          # inferddimproposal::Vector{Float64},
                          dens::Vector{BallTreeDensity},
                          partials::Dict{Int, Vector{BallTreeDensity}};
                          solveKey::Symbol=:default,
                          N::Int=100,
                          dbg::Bool=false  )
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
  inferddimproposal::Vector{Float64}
end

"""
    $SIGNATURES

Calculate the proposals and products on `destvert` using `factors` in factor graph `dfg`.

Notes
- Returns tuple of product and whether full dimensional (=true) or partial (=false).
"""
function predictbelief( dfg::AbstractDFG,
                        destvert::DFGVariable,
                        factors::Vector{<:DFGFactor};
                        N::Int=0,
                        dbg::Bool=false,
                        logger=ConsoleLogger()  )
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

function predictbelief( dfg::AbstractDFG,
                        destvertsym::Symbol,
                        factorsyms::Vector{Symbol};
                        N::Int=0,
                        dbg::Bool=false,
                        logger=ConsoleLogger()  )
  #
  factors = map(fsym -> DFG.getFactor(dfg, fsym), factorsyms)

  vert = DFG.getVariable(dfg, destvertsym)

  # determine the number of particles to draw from the marginal
  nn = N != 0 ? N : size(getVal(vert),2)

  # do the belief prediction
  predictbelief(dfg, vert, factors, N=nn, dbg=dbg, logger=logger)
end

function predictbelief( dfg::AbstractDFG,
                        destvertsym::Symbol,
                        factorsyms::Colon;
                        N::Int=0,
                        dbg::Bool=false,
                        logger=ConsoleLogger() )
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
  pp = AMP.manikde!(pGM, getVariableType(vari) |> getManifolds )

  return pp, dens, partials, lb, sum(inferdim)
end
localProduct(dfg::AbstractDFG, lbl::AbstractString; solveKey::Symbol=:default, N::Int=100, dbg::Bool=false) = localProduct(dfg, Symbol(lbl), solveKey=solveKey, N=N, dbg=dbg)




"""
    $SIGNATURES

Basic wrapper to take local product and then set the value of `sym` in `dfg`.

Notes
- returns `::Tuple{BallTreeDensity, Float64, Vector{Symbol}}`

DevNotes:
- Unknown issue first occurred here near IIF v0.8.4 tag, recorded case at 2020-01-17T15:26:17.673
"""
function localProductAndUpdate!(dfg::AbstractDFG,
                                sym::Symbol,
                                setkde::Bool=true,
                                logger=ConsoleLogger() )
  #
  # calculate new points for sym using existing structure around sym in dfg
  newPts, dens, parts, lbl, infdim = localProduct(dfg, sym, N=getSolverParams(dfg).N, logger=logger)
  # maybe update dfg sym with newly calculated points
  setkde && 0 < size(getPoints(newPts),2) ? setValKDE!(dfg, sym, newPts, false, infdim) : nothing

  return newPts, infdim, lbl
end



#