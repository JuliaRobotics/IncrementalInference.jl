
export calcFactorResidual

function approxConvBelief(
  dfg::AbstractDFG,
  fc::FactorCompute,
  target::Symbol,
  measurement::AbstractVector = Tuple[];
  solveKey::Symbol = :default,
  N::Int = length(measurement),
  nullSurplus::Real = 0,
  skipSolve::Bool = false,
  keepCalcFactor::Union{Nothing, <:Channel} = nothing,
)
  #
  v_trg = getVariable(dfg, target)
  N_ = if N != 0
    N 
  elseif hasState(v_trg, solveKey)
    getNumPts(v_trg; solveKey)
  else
    getSolverParams(dfg).N
  end
  # N = N == 0 ? getNumPts(v_trg; solveKey) : N

  # NOTE approxConv results happen in  duplicate memory destination,  ccw.varValsAll always points directly to variable.VND.val
  pts, ipc = evalFactor(
    dfg, 
    fc, 
    v_trg.label, 
    measurement; 
    solveKey, 
    N = N_, 
    skipSolve, 
    nullSurplus,
    keepCalcFactor
  )

  ## FIXME, bad way to find partial info!!!!
  # Not sufficient to use only observability to determine partial, but is necessary
  # original need is if observability on some coords are zero after a convolution
  len = length(ipc)
  mask = 1e-14 .< abs.(ipc)
  partl = collect(1:len)[mask] 
  
  # is the convolution infoPerCoord full or partial
  statekind = getStateKind(v_trg)
  # FIXME, this if induces type instability via partial
  res = if sum(mask) == getDimension(v_trg)
    # not partial
    HomotopyDensity_legacy(statekind, pts; partial = nothing)
  else
    # is partial
    HomotopyDensity_legacy(statekind, pts; partial = partl)
  end
    
  return res
end

approxConv(w...; kw...) = getPoints(approxConvBelief(w...; kw...), false)

"""
    $SIGNATURES

Calculate the sequential series of convolutions in order as listed by `fctLabels`, and starting from the 
value already contained in the first variable.  

Notes
- `target` must be a variable.
- The ultimate `target` variable must be given to allow path discovery through n-ary factors.
- Fresh starting point will be used if first element in `fctLabels` is a unary `<:AbstractPriorObservation`.
- This function will not change any values in `dfg`, and might have slightly less speed performance to meet this requirement.
- pass in `tfg` to get a recoverable result of all convolutions in the chain.

DevNotes
- TODO strong requirement that this function is super efficient on single factor/variable case!
- FIXME must consolidate with `accumulateFactorMeans`
- TODO `solveKey` not fully wired up everywhere yet
  - tfg gets all the solveKeys inside the source `dfg` variables
  - Consolidate with [`accumulateFactorMeans`](@ref), `approxConvBinary`

Related

[`approxDeconv`](@ref), `findShortestPathDijkstra`
"""
function approxConvBelief(
  dfg::AbstractDFG,
  from::Symbol,
  target::Symbol,
  measurement::AbstractVector = Tuple[];
  solveKey::Symbol = :default,
  N::Int = length(measurement),
  tfg::AbstractDFG = LocalDFG(;solverParams=getSolverParams(dfg)),
  path::AbstractVector{Symbol} = Symbol[],
  skipSolve::Bool = false,
  nullSurplus::Real = 0,
  keepCalcFactor::Union{Nothing, <:Channel} = nothing,
)
  #
  # @assert isVariable(dfg, target) "approxConv(dfg, from, target,...) where `target`=$target must be a variable in `dfg`"

  if from in ls(dfg, target)
    # direct request
    # TODO avoid this allocation for direct cases ( dfg, :x1x2f1, :x2[/:x1] )
    path = Symbol[from; target]
    varLbls = Symbol[target;]
  else
    # must first discover shortest factor path in dfg
    # TODO DFG only supports LocalDFG.findShortestPathDijkstra at the time of writing (DFG v0.10.9)
    path = 0 == length(path) ? findShortestPathDijkstra(dfg, from, target) : path
    @assert path[1] == from "sanity check failing for shortest path function"

    # list of variables
    fctMsk = isFactor.(dfg, path)
    # which factors in the path
    fctLbls = path[fctMsk]
    # must still add
    varLbls = union(lsf.(dfg, fctLbls)...)
    neMsk = exists.(tfg, varLbls) .|> x -> xor(x, true)
    # put the non-existing variables into the temporary graph `tfg`
    # bring all the solveKeys too
    for v in getVariable.(dfg, varLbls[neMsk])
      addVariable!(tfg, v.label, getStateKind(v))
    end
    # variables adjacent to the shortest path should be initialized from dfg
    setdiff(varLbls, path[xor.(fctMsk, true)]) .|>
    x -> initVariable!(tfg, x, getBelief(dfg, x))
  end

  # find/set the starting point
  idxS = 1
  pts = if varLbls[1] == from
    # starting from a variable
    getBelief(dfg, varLbls[1]) |> getPoints
  else
    # chain would start one later
    idxS += 1
    # get the factor
    fct0 = getFactor(dfg, from)
    # get the Matrix{<:Real} of projected points
    pts1Bel = approxConvBelief(
      dfg,
      fct0,
      path[2],
      measurement;
      solveKey,
      N,
      skipSolve,
      nullSurplus,
      keepCalcFactor,
    )
    if length(path) == 2
      return pts1Bel
    end
    getPoints(pts1Bel)
  end
  # didn't return early so shift focus to using `tfg` more intensely
  # FIXME, since AMP v0.15, cannot just set the points, must set HomotopyDensity
  initVariable!(tfg, varLbls[1], pts)

  # do chain of convolutions
  for idx = idxS:length(path)
    if fctMsk[idx]
      # this is a factor path[idx]
      fct = getFactor(dfg, path[idx])
      addFactor!(tfg, fct)
      ptsBel = approxConvBelief(tfg, fct, path[idx + 1]; solveKey, N, skipSolve, keepCalcFactor)
      initVariable!(tfg, path[idx + 1], ptsBel)
    end
  end

  # return target variable values
  return getBelief(tfg, target)
end

"""
    $(SIGNATURES)

Compute proposal belief on `vertid` through `fct` representing some constraint in factor graph.
Always full dimension variable node -- partial constraints will only influence subset of variable dimensions.
The remaining dimensions will keep pre-existing variable values.

Notes
- fulldim is true when "rank-deficient" -- TODO swap to false (or even float)
"""
function calcProposalBelief(
  dfg::AbstractDFG,
  fct::FactorCompute,
  target::Symbol,
  measurement::AbstractVector = Tuple[];
  N::Int = length(measurement),
  solveKey::Symbol = :default,
  nullSurplus::Real = 0,
  dbg::Bool = false,
  keepCalcFactor::Union{Nothing, <:Channel} = nothing,
)
  #
  # assuming it is properly initialized TODO
  proposal = approxConvBelief(dfg, fct, target, measurement; solveKey, N, nullSurplus, keepCalcFactor)

  # _whatP(::HomotopyDensityLive{H, P}) where {H, P} = P
  # @info "calcProposalBelief" getLabel(fct) target _whatP(proposal)

  # return the proposal belief and inferdim, NOTE likely to be changed
  return proposal
end

# specifically the PartialPriorPassThrough dispatch
function calcProposalBelief(
  dfg::AbstractDFG,
  fct::FactorCompute{<:CommonConvWrapper{<:PartialPriorPassThrough}},
  target::Symbol,
  measurement::AbstractVector = Tuple[];
  N::Int = length(measurement),
  solveKey::Symbol = :default,
  nullSurplus::Real = 0,
  dbg::Bool = false,
)
  #

  # density passed through directly from PartialPriorPassThrough.Z
  fctFnc = getObservation(fct)
  proposal = fctFnc.Z.heatmap.densityFnc

  # in case of partial, place the proposal into larger marginal/partial MKD
  proposal_ = if isPartial(fctFnc)
    # oldbel = getBelief(dfg, target, solveKey)
    varType = getStateKind(dfg, target)
    M = getManifold(varType)
    u0 = getPointIdentity(varType)
    # replace(oldbel, proposal)
    antimarginal(M, u0, proposal, Int[fctFnc.partial...])
  else
    proposal
  end

  # return the proposal belief and inferdim, NOTE likely to be changed
  return proposal_
end

"""
    $SIGNATURES

Compute the proposals of a destination vertex for each of `factors` and place the result
as belief estimates in both `dens` and `partials` respectively.

Notes
- TODO: also return if proposals were "dimension-deficient" (aka ~rank-deficient).
"""
function proposalbeliefs!(
  dfg::AbstractDFG,
  destlbl::Symbol,
  factors::AbstractVector, #{<:FactorCompute},
  dens::AbstractVector{<:ApproxManifoldProducts.HomotopyDensity}, # TODO, convert promote to avoid union-abstract vector
  measurement::AbstractVector = Tuple[];
  solveKey::Symbol = :default,
  N::Int = getSolverParams(dfg).N, #maximum([length(getPoints(getBelief(dfg, destlbl, solveKey))); getSolverParams(dfg).N]),
  # how much nullSurplus should be added, see #1517
  nullSurplusAdd::Real = getSolverParams(dfg).nullSurplusAdd,
  dbg::Bool = false,
)
  #

  # populate the full and partial dim containers
  ipcs = Vector{Vector{Float64}}(undef, length(factors))

  # workaround for IIF #1517, additional entropy for sibling factors to target variable if one has multihypo
  nullSrp = zeros(length(factors))
  if any(isMultihypo.(factors))
    # relative sibling factors get nullSurplus
    for (i, f) in enumerate(factors)
      # don't add additional nullSurplus, since its already being done in ExplicitDiscreteMarg!!!  FIXME refactor to common solution
      if isa(getObservation(f), AbstractRelativeObservation) && !isMultihypo(f)
        nullSrp[i] = nullSurplusAdd
      end
    end
  end

  vardim = getDimension(getVariable(dfg, destlbl))
  # get a proposal belief from each factor connected to destlbl
  for (count, fct) in enumerate(factors)
    # need way to convey partial information
    # determine if evaluation is "dimension-deficient" solvable dimension
    # FIXME, update to infoPerCoord
    fct_ipc = ones(vardim) # getFactorSolvableDim(dfg, fct, destlbl, solveKey)
    # convolve or passthrough to get a new proposal
    propBel_ = calcProposalBelief(
      dfg,
      fct,
      destlbl,
      measurement;
      N,
      dbg,
      solveKey,
      nullSurplus = nullSrp[count],
    )
    # partial density
    obs = DFG.getObservation(fct)
    propBel = if isPartial(obs)
      AMP.marginal(propBel_, Int[obs.partial...])
    else
      propBel_
    end
    push!(dens, propBel)
    ipcs[count] = fct_ipc
  end
  # len = maximum(length.(ipcs))
  ipc = zeros(vardim)
  for _ipc in ipcs
    ipc .+= _ipc
  end

  return ipc
end
# group partial dimension factors by selected dimensions -- i.e. [(1,)], [(1,2),(1,2)], [(2,);(2;)]

# WIP, see `_buildGraphByFactorAndTypes!` where pts are full MKD Beliefs, following #1351 
# Legacy use in RoMEPlotting: plotFactor
# function approxConvBelief(fct::AbstractFactorRelative,
#                           varTypes::Union{<:Tuple,<:AbstractVector{<:InstanceType{T}}}, 
#                           mkds::Union{<:Tuple,<:AbstractVector{<:InstanceType{T}}};
#                           tfg::AbstractDFG=_buildGraphByFactorAndTypes!(fct,)
#                           ) where {T <: StateType}
#   #

# end

#
