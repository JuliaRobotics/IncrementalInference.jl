


"""
    $SIGNATURES

Calculate the proposals and products on `destvert` using `factors` in factor graph `dfg`.

Notes
- Returns tuple of product and whether full dimensional (=true) or partial (=false).
- `N` determines the number of samples to draw from the marginal.
- `dens` can contain mixed full and partial dimension `ManifoldKernelDensity` beliefs
"""
function predictbelief( dfg::AbstractDFG,
                        destvert::DFGVariable,
                        factors::Vector{<:DFGFactor};
                        solveKey::Symbol=:default,
                        dens = Vector{ManifoldKernelDensity}(),
                        N::Int=length(getVal(destvert, solveKey=solveKey)),
                        needFreshMeasurements::Bool=true,
                        dbg::Bool=false,
                        logger=ConsoleLogger()  )
  #
  
  # get proposal beliefs
  destvertlabel = destvert.label
  inferdim = proposalbeliefs!(dfg, destvertlabel, factors, dens, solveKey=solveKey, N=N, dbg=dbg)

  # take the product
  oldpts = getBelief(dfg, destvertlabel, solveKey) |> getPoints
  varType = getVariableType(dfg, destvertlabel)
  pGM = AMP.productbelief(oldpts, getManifold(varType), dens, N, dbg=dbg, logger=logger, asPartial=false )

  return pGM, sum(inferdim)
end

predictbelief(dfg::AbstractDFG,
              destlbl::Symbol,
              fctlbls::AbstractVector{Symbol};
              kw... ) = predictbelief(dfg, getVariable(dfg, destlbl), getFactor.(dfg, fctlbls); kw... )
#

predictbelief(dfg::AbstractDFG,
              destlbl::Symbol,
              ::Colon;
              kw... ) = predictbelief(dfg, destlbl, getNeighbors(dfg, destlbl); kw... )
#


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
  # partials = Dict{Any, Vector{ManifoldKernelDensity}}()
  pGM, sinfd = predictbelief(dfg, sym, lb, solveKey=solveKey, logger=logger, dens=dens )

  # make manifold belief from product
  vari = getVariable(dfg, sym)
  pp = AMP.manikde!(getManifold(getVariableType(vari)), pGM )

  return pp, dens, lb, sinfd
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
  newPts, dens, lbl, infdim = localProduct(dfg, sym, solveKey=solveKey, N=getSolverParams(dfg).N, logger=logger)
  # maybe update dfg sym with newly calculated points
  setkde && 0 < length(getPoints(newPts)) ? setValKDE!(dfg, sym, newPts, false, infdim, solveKey=solveKey) : nothing

  return newPts, infdim, lbl
end



#