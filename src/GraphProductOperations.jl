


"""
    $SIGNATURES

Calculate the proposals and products on `destvert` using `factors` in factor graph `dfg`.

Notes
- Returns tuple of product and whether full dimensional (=true) or partial (=false).
"""
function predictbelief( dfg::AbstractDFG,
                        destvert::DFGVariable,
                        factors::Vector{<:DFGFactor};
                        solveKey::Symbol=:default,
                        needFreshMeasurements::Bool=true,
                        N::Int=0,
                        dbg::Bool=false,
                        logger=ConsoleLogger(),
                        dens = Array{ManifoldKernelDensity,1}(),
                        partials = Dict{Any, Vector{ManifoldKernelDensity}}()  )
  #
  
  # determine number of particles to draw from the marginal
  nn = N != 0 ? N : size(getVal(destvert, solveKey=solveKey),2)
  
  # get proposal beliefs
  destvertlabel = destvert.label
  inferdim = proposalbeliefs!(dfg, destvertlabel, factors, dens, partials, solveKey=solveKey, N=nn, dbg=dbg)

  # take the product
  oldpts = getBelief(dfg, destvertlabel, solveKey) |> getPoints
  varType = getVariableType(dfg, destvertlabel)
  pGM = AMP.productbelief(oldpts, getManifold(varType), dens, partials, nn, dbg=dbg, logger=logger )

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
                        dens = Array{ManifoldKernelDensity,1}(),
                        partials = Dict{Any, Vector{ManifoldKernelDensity}}()  )
  #
  factors = getFactor.(dfg, factorsyms)
  vert = getVariable(dfg, destvertsym)

  # determine the number of particles to draw from the marginal
  nn = N != 0 ? N : length(getVal(vert, solveKey=solveKey))

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
                        dens = Array{ManifoldKernelDensity,1}(),
                        partials = Dict{Any, Vector{ManifoldKernelDensity}}() )
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
  dens = Array{ManifoldKernelDensity,1}()
  partials = Dict{Any, Vector{ManifoldKernelDensity}}()
  pGM, sinfd = predictbelief(dfg, sym, lb, solveKey=solveKey, logger=logger, dens=dens, partials=partials)

  # make manifold belief from product
  vari = getVariable(dfg, sym)
  pp = AMP.manikde!(pGM, getVariableType(vari) |> getManifold )

  return pp, dens, partials, lb, sinfd
end
localProduct(dfg::AbstractDFG, lbl::AbstractString; solveKey::Symbol=:default, N::Int=100, dbg::Bool=false) = localProduct(dfg, Symbol(lbl), solveKey=solveKey, N=N, dbg=dbg)




"""
    $SIGNATURES

Basic wrapper to take local product and then set the value of `sym` in `dfg`.

Notes
- returns `::Tuple{ManifoldKernelDensity, Float64, Vector{Symbol}}`

DevNotes:
- Unknown issue first occurred here near IIF v0.8.4 tag, recorded case at 2020-01-17T15:26:17.673
"""
function localProductAndUpdate!(dfg::AbstractDFG,
                                sym::Symbol,
                                setkde::Bool=true,
                                logger=ConsoleLogger();
                                solveKey::Symbol=:default )
  #
  # calculate new points for sym using existing structure around sym in dfg
  newPts, dens, parts, lbl, infdim = localProduct(dfg, sym, solveKey=solveKey, N=getSolverParams(dfg).N, logger=logger)
  # maybe update dfg sym with newly calculated points
  setkde && 0 < size(getPoints(newPts),2) ? setValKDE!(dfg, sym, newPts, false, infdim, solveKey=solveKey) : nothing

  return newPts, infdim, lbl
end



#