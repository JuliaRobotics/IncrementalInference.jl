

function _partialProducts!(pGM, partials, manis; useExisting::Bool=false)
  if useExisting
    # include previous calcs
    for (dimnum,pp) in partials
      dimv = [dimnum...]
      push!(pp, AMP.manikde!(pGM[dimv,:], (manis[dimv]...,) ))
    end
  end

  # do each partial dimension individually
  for (dimnum,pp) in partials
    dimv = [dimnum...]
    pGM[dimv,:] = AMP.manifoldProduct(pp, (manis[dimv]...,), Niter=1) |> getPoints
  end
end

# """
#     $(SIGNATURES)

# Multiply various full and partial dimension proposal densities.

# DevNotes
# - FIXME consolidate partial and full product AMP API, relates to #1010
# - TODO better consolidate with full dimension product
# - TODO -- reuse memory rather than rand here
# """
# function prodmultiplefullpartials(dens::Vector{BallTreeDensity},
#                                   partials::Dict{Any, Vector{BallTreeDensity}},
#                                   Ndims::Int,
#                                   N::Int,
#                                   manis::Tuple;
#                                   useExisting::Bool=false )
#   #
#   # calculate products over all dimensions, legacy proposals held in `dens` vector
#   pGM = AMP.manifoldProduct(dens, manis, Niter=1) |> getPoints

#   _partialProducts!(pGM, partials, manis, useExisting=useExisting)

#   return pGM
# end



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
                        dens::Vector{<:BallTreeDensity},
                        partials::Dict{Any, <:AbstractVector{<:BallTreeDensity}},
                        N::Int;
                        dbg::Bool=false,
                        logger=ConsoleLogger()  )
  #
  vert = DFG.getVariable(dfg, vertlabel)
  manis = getVariableType(vert) |> getManifolds
  pGM = Array{Float64,2}(undef, 0,0)
  lennonp, lenpart = length(dens), length(partials)
  denspts = getPoints(getBelief(vert))
  Ndims = size(denspts,1)
  with_logger(logger) do
    @info "[$(lennonp)x$(lenpart)p,d$(Ndims),N$(N)],"
  end

  dens_ = 0 < lennonp ? dens : [deepcopy(getBelief(vert));]
    # # resize for #1013
    # if size(denspts,2) < N
    #   pGM = zeros(size(denspts,1),N)
    #   pGM[:,1:size(denspts,2)] .= denspts
    # else
    #   pGM = deepcopy(denspts)
    # end

  if 0 < lennonp # || (lennonp == 0 && 0 < lenpart)
    # calculate products over all dimensions, legacy proposals held in `dens` vector
    pGM = AMP.manifoldProduct(dens, manis, Niter=1) |> getPoints
    # multiple non-partials
    _partialProducts!(pGM, partials, manis, useExisting=true)
    # pGM = prodmultiplefullpartials(dens_, partials, Ndims, N, manis, useExisting=true)
  elseif lennonp == 0 && 0 < lenpart
    # only partials, must get other existing values for vertlabel from dfg
    pGM = deepcopy(denspts)
    _partialProducts!(pGM, partials, manis; useExisting=false)
    #   # do each partial dimension individually
    # for (dimnum,pp) in partials
    #   pGM[dimnum,:] = AMP.manifoldProduct(pp, (manis[dimnum],), Niter=1) |> getPoints
    # end
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

Calculate the proposals and products on `destvert` using `factors` in factor graph `dfg`.

Notes
- Returns tuple of product and whether full dimensional (=true) or partial (=false).
"""
function predictbelief( dfg::AbstractDFG,
                        destvert::DFGVariable,
                        factors::Vector{<:DFGFactor};
                        solveKey::Symbol=:defaul,
                        needFreshMeasurements::Bool=true,
                        N::Int=0,
                        dbg::Bool=false,
                        logger=ConsoleLogger(),
                        dens = Array{BallTreeDensity,1}(),
                        partials = Dict{Any, Vector{BallTreeDensity}}()  )
  #
  
  # determine number of particles to draw from the marginal
  nn = N != 0 ? N : size(getVal(destvert),2)
  
  # get proposal beliefs
  destvertlabel = destvert.label
  inferdim = proposalbeliefs!(dfg, destvertlabel, factors, dens, partials, N=nn, dbg=dbg)

  # take the product
  pGM = productbelief(dfg, destvertlabel, dens, partials, nn, dbg=dbg, logger=logger )

  return pGM, sum(inferdim)
end

function predictbelief( dfg::AbstractDFG,
                        destvertsym::Symbol,
                        factorsyms::AbstractVector{Symbol};
                        solveKey::Symbol=:default,
                        needFreshMeasurements::Bool=true,
                        N::Int=0,
                        dbg::Bool=false,
                        logger=ConsoleLogger(),
                        dens = Array{BallTreeDensity,1}(),
                        partials = Dict{Any, Vector{BallTreeDensity}}()  )
  #
  factors = getFactor.(dfg, factorsyms)
  vert = getVariable(dfg, destvertsym)

  # determine the number of particles to draw from the marginal
  nn = N != 0 ? N : size(getVal(vert, solveKey=solveKey),2)

  # do the belief prediction
  predictbelief(dfg, vert, factors, solveKey=solveKey, needFreshMeasurements=needFreshMeasurements, N=nn, dbg=dbg, logger=logger, dens=dens, partials=partials)
end

function predictbelief( dfg::AbstractDFG,
                        destvertsym::Symbol,
                        ::Colon;
                        solveKey::Symbol=:default,
                        needFreshMeasurements::Bool=true,
                        N::Int=0,
                        dbg::Bool=false,
                        logger=ConsoleLogger(),
                        dens = Array{BallTreeDensity,1}(),
                        partials = Dict{Any, Vector{BallTreeDensity}}() )
  #
  predictbelief(dfg, destvertsym, getNeighbors(dfg, destvertsym), solveKey=solveKey, needFreshMeasurements=needFreshMeasurements, N=N, dbg=dbg, logger=logger, dens=dens, partials=partials )
end

"""
    $(SIGNATURES)

Using factor graph object `dfg`, project belief through connected factors
(convolution with likelihood) to variable `sym` followed by a approximate functional product.

Return: product belief, full proposals, partial dimension proposals, labels
"""
function localProduct(dfg::AbstractDFG,
                      sym::Symbol;
                      solveKey::Symbol=:default,
                      N::Int=100,
                      dbg::Bool=false,
                      logger=ConsoleLogger() )
  #
  # vector of all neighbors as Symbols
  lb = getNeighbors(dfg, sym)

  # # get proposal beliefs
  dens = Array{BallTreeDensity,1}()
  partials = Dict{Any, Vector{BallTreeDensity}}()
  pGM, sinfd = predictbelief(dfg, sym, lb, solveKey=solveKey, logger=logger, dens=dens, partials=partials)

  # make manifold belief from product
  vari = getVariable(dfg, sym)
  pp = AMP.manikde!(pGM, getVariableType(vari) |> getManifolds )

  return pp, dens, partials, lb, sinfd
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