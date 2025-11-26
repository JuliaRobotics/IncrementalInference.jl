# Factor Graph OS type utilities
#  IIF methods should direclty detect extended types from user import
# of convert in their namespace

# FIXME, upgrade to AMP instead
KDE.getPoints(dfg::AbstractDFG, lbl::Symbol) = getBelief(dfg, lbl) |> getPoints

clampStringLength(st::AbstractString, len::Int = 5) = st[1:minimum([len; length(st)])]

function clampBufferString(
  st::AbstractString,
  max::Int,
  len::Int = minimum([max, length(st)]),
)
  @assert 0 <= max "max must be greater or equal to zero"
  st = clampStringLength(st, len)
  for i = len:(max - 1)
    st *= " "
  end
  return st
end

"""
    $SIGNATURES

Extract contiguous string of numbers at end of a label`::Symbol` -- e.g. `:x45_4` --> "4".  
Returns `(string, suffix_substring)`

Related

[`incrSuffix`](@ref)
"""
function _getSuffix(lbl::Symbol; pattern::Regex = r"\d+")
  slbl = string(lbl)
  phrase_ = slbl |> reverse |> x -> match(pattern, x).match
  return slbl, reverse(phrase_)
end

"""
    $SIGNATURES

Utility for incrementing or decrementing suffix numbers in DFG variable labels, e.g.
```julia
incrSuffix(:x45_4)
# returns :x45_5

incrSuffix(:x45, +3)
# returns :x48

incrSuffix(:x45_4, -1)
# returns :x45_3
```

Notes
- Change `pattern::Regex=r"\\d+"` for alternative behaviour.
"""
function incrSuffix(lbl::Symbol, val::Integer = +1; pattern::Regex = r"\d+")
  slbl, phrase = _getSuffix(lbl; pattern = pattern)
  nint = phrase |> x -> (parse(Int, x) + val)
  prefix = slbl[1:(end - length(phrase))]
  return Symbol(prefix, nint)
end


_getCCW(fct::FactorCompute) = DFG.getCache(fct) #getState(fct) |> _getCCW
_getCCW(dfg::AbstractDFG, lbl::Symbol) = DFG.getCache(getFactor(dfg, lbl)) #getFactor(dfg, lbl) |> _getCCW

DFG.getObservation(ccw::CommonConvWrapper) = ccw.usrfnc!

_getZDim(ccw::CommonConvWrapper) = getManifold(ccw) |> manifold_dimension # ccw.zDim
# TODO is MsgPrior piggy backing zdim on inferdim???
_getZDim(ccw::CommonConvWrapper{<:MsgPrior}) = length(ccw.usrfnc!.infoPerCoord) # ccw.usrfnc!.inferdim

_getZDim(fct::FactorCompute) = _getCCW(fct) |> _getZDim

DFG.getDimension(fct::FactorCompute) = _getZDim(fct)

"""
    $SIGNATURES

Return the manifold on which this ManifoldKernelDensity is defined.

DevNotes
- TODO currently ignores the .partial aspect (captured in parameter `L`)
"""
function getManifold(
  mkd::ManifoldKernelDensity{M, B, Nothing},
  asPartial::Bool = false,
) where {M, B}
  return mkd.manifold
end
function getManifold(
  mkd::ManifoldKernelDensity{M, B, L},
  asPartial::Bool = false,
) where {M, B, L <: AbstractVector}
  return asPartial ? mkd.manifold : getManifoldPartial(mkd.manifold, mkd._partial)
end

"""
    $TYPEDSIGNATURES

Return the number of dimensions this factor vertex `fc` influences.

DevNotes
- TODO document how this function handles partial dimensions
  - Currently a factor manifold is just what the measurement provides (i.e. bearing only would be dimension 1)
"""
getFactorDim(w...) = getDimension(w...)
getFactorDim(fg::AbstractDFG, fctid::Symbol) = getFactorDim(getFactor(fg, fctid))

# extend convenience function (Matrix or Vector{P})
function manikde!(
  variableType::Union{InstanceType{<:StateType}, InstanceType{<:AbstractObservation}},
  pts::AbstractVector{P};
  kw...,
) where {P <: Union{<:AbstractArray, <:Number, <: ArrayPartition}}
  #
  M = getManifold(variableType)
  # @info "pts" P typeof(pts[1]) pts[1]
  infoPerCoord = ones(AMP.getNumberCoords(M, pts[1]))
  return AMP.manikde!(M, pts; infoPerCoord, kw...)
end

function manikde!(
  varT::InstanceType{<:StateType},
  pts::AbstractVector{<:Tuple};
  kw...,
)
  #
  return manikde!(varT, (t -> ArrayPartition(t...)).(pts); kw...)
end

"""
    $SIGNATURES

Return params.N measurement samples for a factor in `<:AbstractDFG`.
"""
function getMeasurements(dfg::AbstractDFG, fsym::Symbol, N::Int = getSolverParams(dfg).N)
  return sampleFactor(dfg, fsym, N)
end

"""
    $SIGNATURES

Get the folder location where debug and solver information is recorded for a particular factor graph.
"""
getLogPath(opt::SolverParams) = opt.logpath
getLogPath(dfg::AbstractDFG) = getSolverParams(dfg) |> getLogPath

"""
    $SIGNATURES

Append `str` onto factor graph log path as convenience function.
"""
joinLogPath(opt::SolverParams, str...) = joinpath(getLogPath(opt), str...)
joinLogPath(dfg::AbstractDFG, str...) = joinLogPath(getSolverParams(dfg), str...)

"""
    $(SIGNATURES)

Set variable(s) `sym` of factor graph to be marginalized -- i.e. not be updated by inference computation.
"""
function setfreeze!(dfg::AbstractDFG, sym::Symbol, solveKey::Symbol = :default)
  if !isInitialized(dfg, sym)
    @warn "Vertex $(sym) is not initialized, and won't be frozen at this time."
    return nothing
  end
  vert = DFG.getVariable(dfg, sym)
  data = getState(vert, solveKey)
  data.marginalized = true
  return nothing
end
function setfreeze!(dfg::AbstractDFG, syms::Vector{Symbol}, solveKey::Symbol = :default)
  for sym in syms
    setfreeze!(dfg, sym, solveKey)
  end
end

"""
    $(SIGNATURES)

Freeze nodes that are older than the quasi fixed-lag length defined by `fg.qfl`, according to `fg.fifo` ordering.

Future:
- Allow different freezing strategies beyond fifo.
"""
function fifoFreeze!(dfg::AbstractDFG)
  if DFG.getSolverParams(dfg).qfl == 0
    @warn "Quasi fixed-lag is enabled but QFL horizon is zero. Please set a valid window with FactoGraph.qfl"
  end

  # the fifo history
  tofreeze = DFG.getAddHistory(dfg)[1:(end - DFG.getSolverParams(dfg).qfl)]

  # check that the variable to freeze exists fix issue #966
  filter!(v -> exists(dfg, v), tofreeze)

  if length(tofreeze) == 0
    @info "[fifoFreeze] QFL - no nodes to freeze."
    return nothing
  end
  @info "[fifoFreeze] QFL - Freezing nodes $(tofreeze[1]) -> $(tofreeze[end])."
  setfreeze!(dfg, tofreeze)
  return nothing
end

DFG.getPoint(typ::StateType, w...; kw...) = getPoint(typeof(typ), w...; kw...)
function DFG.getCoordinates(typ::StateType, w...; kw...)
  return getCoordinates(typeof(typ), w...; kw...)
end

# WIP
# _getMeasurementRepresentation(::AbstractPriorObservation, coord::AbstractVector{<:Number}) = 

"""
    $SIGNATURES

Get the ParametricPointEstimates---based on full marginal belief estimates---of a variable in the distributed factor graph.
Calculate new Parametric Point Estimates for a given variable.


DevNotes
- TODO update for manifold subgroups.
- TODO standardize after AMP3D
"""
function calcMeanMaxSuggested(
  vari::VariableCompute,
  solveKey::Symbol = :default
)
  varType = getStateKind(vari)
  P = getBelief(vari, solveKey)
  maniDef = convert(MB.AbstractManifold, varType)
  manis = AMP._manifoldtuple(maniDef) # LEGACY, TODO REMOVE
  ops = buildHybridManifoldCallbacks(manis)
  Pme = calcMean(P)  # getKDEMean(P) #, addop=ops[1], diffop=ops[2]

  # returns coordinates at identify
  Pma = getKDEMax(P; addop = ops[1], diffop = ops[2])
  # calculate point

  ## TODO use getCoordinates for now (IIF v0.25)
  Pme_ = getCoordinates(varType, Pme)
  # Pma_ = getCoordinates(M,Pme)
 
  return (
    mean=Pme_, 
    max=Pma, 
    suggested=Pme_, 
  )
end

function calcMeanMaxSuggested(
  dfg::AbstractDFG,
  label::Symbol,
  solveKey::Symbol = :default,
)
  return calcMeanMaxSuggested(getVariable(dfg, label), solveKey)
end

"""
    $SIGNATURES

Return bool on whether a certain factor has user defined multihypothesis.

Related

[`getMultihypoDistribution`](@ref)
"""
isMultihypo(fct::FactorCompute) = isa(_getCCW(fct).hyporecipe.hypotheses, Distribution)

"""
    $SIGNATURES

Return the categorical distributed used for multihypothesis selection in a factor.

Related

isMultihypo
"""
getMultihypoDistribution(fct::FactorCompute) = _getCCW(fct).hyporecipe.hypotheses

"""
    $SIGNATURES

Free all variables from marginalization.
"""
function dontMarginalizeVariablesAll!(fgl::AbstractDFG)
  fgl.solverParams.isfixedlag = false
  fgl.solverParams.qfl = (2^(Sys.WORD_SIZE - 1) - 1)
  fgl.solverParams.limitfixeddown = false
  for sym in ls(fgl)
    setMarginalized!(fgl, sym, false)
  end
  return nothing
end

"""
    $SIGNATURES

Free all variables from marginalization.

Related

dontMarginalizeVariablesAll!
"""
function unfreezeVariablesAll!(fgl::AbstractDFG)
  return dontMarginalizeVariablesAll!(fgl)
end

# WIP
# function resetSolvableAllExcept!(dfg::AbstractDFG,
#                                   fltr::NothingUnion{Regex}=nothing)
#   #
#   unfreezeVariablesAll!(dfg)
# end

"""
    $SIGNATURES

Reset initialization flag on all variables in `::AbstractDFG`.

Notes
- Numerical values remain, but inference will overwrite since init flags are now `false`.
"""
function resetVariableAllInitializations!(fgl::AbstractDFG, solveKey::Symbol = :default)
  vsyms = ls(fgl)
  for sym in vsyms
    setVariableInitialized!(getVariable(fgl, sym), solveKey, false)
  end
  return nothing
end

"""
    $SIGNATURES

Enable defaults for fixed-lag-like operation by using smart message passing on the tree.

Notes:
- These are only default settings, and can be modified in each use case scenario.
- Default does not update downsolve through to leaves of the tree.
"""
function defaultFixedLagOnTree!(
  dfg::AbstractDFG,
  len::Int = 30;
  limitfixeddown::Bool = true,
)
  #
  getSolverParams(dfg).isfixedlag = true
  getSolverParams(dfg).qfl = len
  getSolverParams(dfg).limitfixeddown = limitfixeddown
  return getSolverParams(dfg)
end

"""
    $SIGNATURES

Find all factors that go `from` variable to any other complete variable set within `between`.

Notes
- Developed for downsolve in CSM, expanding the cliqSubFg to include all frontal factors.
"""
function findFactorsBetweenFrom(
  dfg::G,
  between::Vector{Symbol},
  from::Symbol,
) where {G <: AbstractDFG}
  # get all associated factors
  allfcts = ls(dfg, from)

  # remove candidates with neighbors outside between with mask
  mask = ones(Bool, length(allfcts))
  i = 0
  for fct in allfcts
    i += 1
    # check if immediate neighbors are all in the `between` list
    immnei = ls(dfg, fct)
    if length(immnei) != length(intersect(immnei, between))
      mask[i] = false
    end
  end

  # return only masked factors
  return allfcts[mask]
end

"""
    $SIGNATURES

Return list of factors which depend only on variables in variable list in factor
graph -- i.e. among variables.

Notes
-----
* `unused::Bool=true` will disregard factors already used -- i.e. disregard where `potentialused=true`
"""
function getFactorsAmongVariablesOnly(
  dfg::G,
  varlist::Vector{Symbol};
  unused::Bool = true,
) where {G <: AbstractDFG}
  # collect all factors attached to variables
  prefcts = Symbol[]
  for var in varlist
    union!(prefcts, DFG.ls(dfg, var))
  end

  almostfcts = Symbol[]
  if unused
    # now check if those factors have already been added
    for fct in prefcts
      vert = DFG.getFactor(dfg, fct)
      if !vert.state.potentialused
        push!(almostfcts, fct)
      end
    end
  else
    almostfcts = prefcts
  end

  # Select factors that have all variables in this clique var list
  usefcts = Symbol[]
  for fct in almostfcts
    if length(setdiff(listNeighbors(dfg, fct), varlist)) == 0
      push!(usefcts, fct)
    end
  end

  return usefcts
end

## ============================================================================
# Starting integration with Manifolds.jl, via ApproxManifoldProducts.jl first
## ============================================================================

"""
    $SIGNATURES
Fetch and unpack JSON dictionary stored as a data blob.
"""
function fetchDataJSON(dfg::AbstractDFG, varsym::Symbol, lbl::Symbol)
  gde, rawData = getData(dfg, varsym, lbl)
  if gde.mimeType == "application/json/octet-stream"
    JSON3.read(IOBuffer(rawData))
  else
    error("Unknown JSON Blob format $(gde.mimeType)")
  end
end

#
