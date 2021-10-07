

"""
    $SIGNATURES

Calculate the proposals and products on `destvert` using `factors` in factor graph `dfg`.

Notes
- Returns tuple of product and whether full dimensional (=true) or partial (=false).
- `N` determines the number of samples to draw from the marginal.
- `dens` can contain mixed full and partial dimension `ManifoldKernelDensity` beliefs

Related

[`approxConv`](@ref), [`proposalbeliefs!`](@ref), [`AMP.manifoldProduct`](@ref)
"""
function propagateBelief( dfg::AbstractDFG,
                          destvar::DFGVariable,
                          factors::AbstractVector{<:DFGFactor};
                          solveKey::Symbol=:default,
                          dens = Vector{ManifoldKernelDensity}(),
                          N::Int=getSolverParams(dfg).N, #maximum([length(getPoints(getBelief(destvar, solveKey))); getSolverParams(dfg).N]),
                          needFreshMeasurements::Bool=true,
                          dbg::Bool=false,
                          logger=ConsoleLogger()  )
  #
  
  # get proposal beliefs
  destlbl = getLabel(destvar)
  ipc = proposalbeliefs!(dfg, destlbl, factors, dens, solveKey=solveKey, N=N, dbg=dbg)
  
  # @show dens[1].manifold

  # make sure oldpts has right number of points
  oldBel = getBelief(dfg, destlbl, solveKey)
  oldpts = if Npts(oldBel) == N
    getPoints(oldBel)
  else
    sample(oldBel, N)[1]
  end
  
  # few more data requirements
  varType = getVariableType(destvar)
  M = getManifold(varType)
  # @info "BUILDING MKD" varType M isPartial.(dens)
  
  # take the product
  mkd = AMP.manifoldProduct(dens, M, Niter=1, oldPoints=oldpts, N=N)
  
  # @info "GOT" mkd.manifold
  
  return mkd, ipc
end

"""
    $SIGNATURES

This is an old function name that will be replaced by [`propagateBelief`](@ref).
"""
function predictbelief( dfg::AbstractDFG,
                        destvert::DFGVariable,
                        factors::AbstractVector{<:DFGFactor};
                        asPartial::Bool=false,
                        kw... )
  #
  # new
  mkd, ifd = propagateBelief(dfg,destvert,factors;kw...)

  # legacy interface
  return getPoints(mkd, asPartial), ifd
end

predictbelief(dfg::AbstractDFG,
              destlbl::Symbol,
              fctlbls::AbstractVector{Symbol};
              kw... ) = predictbelief(dfg, getVariable(dfg, destlbl), map(x->getFactor(dfg, x), fctlbls); kw... )
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
                      N::Int=getSolverParams(dfg).N, #maximum([length(getPoints(getBelief(dfg, sym, solveKey))); getSolverParams(dfg).N]),
                      dbg::Bool=false,
                      logger=ConsoleLogger() )
  #
  # vector of all neighbors as Symbols
  lb = getNeighbors(dfg, sym)

  # store proposal beliefs, TODO replace Abstract with concrete type
  dens = Vector{ManifoldKernelDensity}()
  
  mkd, sinfd = propagateBelief(dfg, getVariable(dfg, sym), map(x->getFactor(dfg, x), lb); solveKey=solveKey, logger=logger, dens=dens, N=N )
  
  return mkd, dens, lb, sinfd
end
localProduct(dfg::AbstractDFG, lbl::AbstractString; kw...) = localProduct(dfg, Symbol(lbl); kw...)




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
  newPts, dens, lbl, ipc = localProduct(dfg, sym, solveKey=solveKey, N=getSolverParams(dfg).N, logger=logger)
  # maybe update dfg sym with newly calculated points
  setkde && 0 < length(getPoints(newPts)) ? setValKDE!(dfg, sym, newPts, false, ipc, solveKey=solveKey) : nothing

  return newPts, ipc, lbl
end



#