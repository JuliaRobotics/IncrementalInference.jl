
"""
    $SIGNATURES

Calculate the proposals and products on `destvert` using `factors` in factor graph `dfg`.

Notes
- Returns tuple of product and whether full dimensional (=true) or partial (=false).
- `N` determines the number of samples to draw from the marginal.
- `dens` can contain mixed full and partial dimension `ManifoldKernelDensity` beliefs

Related

[`approxConvBelief`](@ref), [`proposalbeliefs!`](@ref), [`AMP.manifoldProduct`](@ref)
"""
function propagateBelief(
  dfg::AbstractDFG,
  destvar::VariableCompute,
  factors::AbstractVector; #{<:FactorCompute};
  solveKey::Symbol = :default,
  dens::AbstractVector{<:ApproxManifoldProducts.HomotopyDensity} = Vector{HomotopyDensityLive}(), # TODO, abstract requires dynamic dispatch (slow)
  N::Integer = getSolverParams(dfg).N,
  needFreshMeasurements::Bool = true,
  dbg::Bool = false,
  logger = ConsoleLogger(),
  asPartial::Bool=false,
)
  # get proposal beliefs
  destlbl = getLabel(destvar)
  _observability = proposalbeliefs!(dfg, destlbl, factors, dens; solveKey, N, dbg)

  # # make sure oldPoints vector has right length
  #   # oldBel = getBelief(dfg, destlbl, solveKey; newbw = false)
  #   # _pts = getPoints(oldBel, false)
  # oldpts = DistributedFactorGraphs.refPoints(DistributedFactorGraphs.getState(destvar, solveKey))
  # if N != length(oldpts)
  #   resize!(oldpts, N)
  # end

  # # few more data requirements
  # varType = getStateKind(destvar)
  
  # take the product
  hode = manifoldProduct(
    dens;
    MC = 1,
    N,
  )

  _whatP(::HomotopyDensityLive{H, P}) where {H, P} = P
  # @info "propagateBelief" getLabel(destvar) _whatP(hode) string(_whatP.(dens))

  return hode, _observability
end

function propagateBelief(
  dfg::AbstractDFG,
  destlbl::Symbol,
  fctlbls::AbstractVector{Symbol};
  kw...,
)
  return propagateBelief(
    dfg,
    getVariable(dfg, destlbl),
    map(x -> getFactor(dfg, x), fctlbls);
    kw...,
  )
end
#

propagateBelief(dfg::AbstractDFG, destlbl::Symbol, ::Colon; kw...) = propagateBelief(dfg, destlbl, listNeighbors(dfg, destlbl); kw...)



"""
    $(SIGNATURES)

Using factor graph object `dfg`, project belief through connected factors
(convolution with likelihood) to variable `sym` followed by a approximate functional product.

Return: product belief, full proposals, partial dimension proposals, labels
"""
function localProduct(
  dfg::AbstractDFG,
  sym::Symbol;
  solveKey::Symbol = :default,
  N::Int = getSolverParams(dfg).N, #maximum([length(getPoints(getBelief(dfg, sym, solveKey))); getSolverParams(dfg).N]),
  dbg::Bool = false,
  logger = ConsoleLogger(),
)
  #
  # vector of all neighbors as Symbols
  lb = listNeighbors(dfg, sym)

  # store proposal beliefs, TODO replace Abstract with concrete type
  dens = Vector{HomotopyDensity}()

  fcts = map(x -> getFactor(dfg, x), lb)
  mkd, sinfd = propagateBelief(
    dfg,
    getVariable(dfg, sym),
    fcts;
    solveKey = solveKey,
    logger = logger,
    dens = dens,
    N = N,
  )

  return mkd, dens, lb, sinfd
end
function localProduct(dfg::AbstractDFG, lbl::AbstractString; kw...)
  return localProduct(dfg, Symbol(lbl); kw...)
end

"""
    $SIGNATURES

Basic wrapper to take local product and then set the value of `sym` in `dfg`.

Notes
- returns `::Tuple{ManifoldKernelDensity, Float64, Vector{Symbol}}`

DevNotes:
- Unknown issue first occurred here near IIF v0.8.4 tag, recorded case at 2020-01-17T15:26:17.673
"""
function localProductAndUpdate!(
  dfg::AbstractDFG,
  sym::Symbol,
  setkde::Bool = true,
  logger = ConsoleLogger();
  solveKey::Symbol = :default,
)
  #
  # calculate new points for sym using existing structure around sym in dfg
  newPts, dens, lbl, ipc =
    localProduct(dfg, sym; solveKey = solveKey, N = getSolverParams(dfg).N, logger = logger)
  # maybe update dfg sym with newly calculated points
  if setkde && 0 < length(getPoints(newPts))
    setValKDE!(dfg, sym, newPts, false, ipc; solveKey = solveKey)
  else
    nothing
  end

  return newPts, ipc, lbl
end

#
