
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
  destvar::DFGVariable,
  factors::AbstractVector; #{<:DFGFactor};
  solveKey::Symbol = :default,
  dens::AbstractVector{<:ManifoldKernelDensity} = Vector{ManifoldKernelDensity}(), # TODO, abstract requires dynamic dispatch (slow)
  N::Integer = getSolverParams(dfg).N,
  needFreshMeasurements::Bool = true,
  dbg::Bool = false,
  logger = ConsoleLogger(),
  asPartial::Bool=false,
)
  #

  # get proposal beliefs
  destlbl = getLabel(destvar)
  ipc = proposalbeliefs!(dfg, destlbl, factors, dens; solveKey, N, dbg)

  # @show dens[1].manifold

  # make sure oldPoints vector has right length
  oldBel = getBelief(dfg, destlbl, solveKey)
  _pts = getPoints(oldBel, false)
  oldPoints = if Npts(oldBel) < N
    nn = N - length(_pts) # should be larger than 0
    _pts_, = sample(oldBel, nn)
    vcat(_pts, _pts_)
  else
    _pts[1:N]
  end

  # few more data requirements
  varType = getVariableType(destvar)
  M = getManifold(varType)
  # @info "BUILDING MKD" varType M isPartial.(dens)
  
  # take the product
  mkd = AMP.manifoldProduct(
    dens,
    M;
    Niter = 1,
    oldPoints,
    N,
    u0 = getPointDefault(varType),
  )

  # @info "GOT" mkd.manifold
  return mkd, ipc
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
  dens = Vector{ManifoldKernelDensity}()

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
