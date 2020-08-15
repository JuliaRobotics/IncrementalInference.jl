# Factor Graph OS type utilities
#  IIF methods should direclty detect extended types from user import
# of convert in their namespace

import DistributedFactorGraphs: AbstractPointParametricEst, loadDFG


export getPPESuggestedAll, findVariablesNear, defaultFixedLagOnTree!
export loadDFG


# export setSolvable!

manikde!(pts::AbstractArray{Float64,2}, vartype::Union{InstanceType{InferenceVariable}, InstanceType{FunctorInferenceType}}) = manikde!(pts, getManifolds(vartype))
manikde!(pts::AbstractArray{Float64,1}, vartype::Type{ContinuousScalar}) = manikde!(reshape(pts,1,:), getManifolds(vartype))

# extend convenience function
function manikde!(pts::AbstractArray{Float64,2},
  bws::Vector{Float64},
  softtype::Union{InstanceType{InferenceVariable}, InstanceType{FunctorInferenceType}}  )
#
manikde!(pts, bws, getManifolds(softtype))
end


"""
    $SIGNATURES

Return N=100 measurement samples for a factor in `<:AbstractDFG`.
"""
function getMeasurements(dfg::AbstractDFG, fsym::Symbol, N::Int=100)
  fnc = getFactorFunction(dfg, fsym)
  # getSample(fnc, N)
  Xi = (v->getVariable(dfg, v)).(getVariableOrder(dfg, fsym))
  freshSamples(fnc, N)
end

"""
    $SIGNATURES

Get graph node (variable or factor) dimension.
"""
getDimension(vartype::InferenceVariable) = vartype.dims #TODO Deprecate
getDimension(vartype::Type{<:InferenceVariable}) = getDimension(vartype())
getDimension(var::DFGVariable) = getDimension(getSofttype(var))
getDimension(fct::DFGFactor) = getSolverData(fct).fnc.zDim

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
function setfreeze!(dfg::AbstractDFG, sym::Symbol)
  if !isInitialized(dfg, sym)
    @warn "Vertex $(sym) is not initialized, and won't be frozen at this time."
    return nothing
  end
  vert = DFG.getVariable(dfg, sym)
  data = getSolverData(vert)
  data.ismargin = true
  nothing
end
function setfreeze!(dfg::AbstractDFG, syms::Vector{Symbol})
  for sym in syms
    setfreeze!(dfg, sym)
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
  tofreeze = DFG.getAddHistory(dfg)[1:(end-DFG.getSolverParams(dfg).qfl)]
  if length(tofreeze) == 0
      @info "[fifoFreeze] QFL - no nodes to freeze."
      return nothing
  end
  @info "[fifoFreeze] QFL - Freezing nodes $(tofreeze[1]) -> $(tofreeze[end])."
  setfreeze!(dfg, tofreeze)
  nothing
end



"""
    $SIGNATURES

Get the ParametricPointEstimates---based on full marginal belief estimates---of a variable in the distributed factor graph.

DevNotes
- TODO update for manifold subgroups.

Related

getVariablePPE, setVariablePosteriorEstimates!, getVariablePPE!
"""
function calcVariablePPE(var::DFGVariable,
                         softt::InferenceVariable;
                         solveKey::Symbol=:default,
                         method::Type{MeanMaxPPE}=MeanMaxPPE  )
  #
  P = getKDE(var, solveKey)
  manis = getManifolds(softt) # getManifolds(vnd)
  ops = buildHybridManifoldCallbacks(manis)
  Pme = getKDEMean(P) #, addop=ops[1], diffop=ops[2]
  Pma = getKDEMax(P, addop=ops[1], diffop=ops[2])
  suggested = zeros(getDimension(var))
  # TODO standardize after AMP3D
  @assert length(manis) == getDimension(var)
  for i in 1:length(manis)
    mani = manis[i]
    if mani == :Euclid
      suggested[i] = Pme[i]
    elseif mani == :Circular
      suggested[i] = Pma[i]
    else
      error("Unknown manifold to find PPE, $softt, $mani")
    end
  end
  MeanMaxPPE(solveKey, suggested, Pma, Pme, now())
end


calcVariablePPE(var::DFGVariable; method::Type{<:AbstractPointParametricEst}=MeanMaxPPE, solveKey::Symbol=:default) = calcVariablePPE(var, getSofttype(var), method=method, solveKey=solveKey)

function calcVariablePPE(dfg::AbstractDFG,
                         sym::Symbol;
                         method::Type{<:AbstractPointParametricEst}=MeanMaxPPE,
                         solveKey::Symbol=:default )
  #
  var = getVariable(dfg, sym)
  calcVariablePPE(var, getSofttype(var), method=method, solveKey=solveKey)
end


"""
    $SIGNATURES

Return interger index of desired variable element.

Example
-------
```julia
pp = RoME.Point2()
getIdx(pp, :posY) # = 2
```

Internal Notes
--------------
- uses number i < 100 for index number, and
- uses +100 offsets to track the minibatch number of the requested dimension
"""
function getIdx(pp::Tuple,
                sym::Symbol,
                i::Int=0)
  #
  i-=100
  for p in pp
    i,j = getIdx(p, sym, i)
    if i > 0
      return i, j
    end
  end
  return i,-1
end
getIdx(pp::Symbol, sym::Symbol, i::Int=0) = pp==sym ? (abs(i)%100+1, div(abs(i)-100,100)) : (i-1, div(abs(i)-100,100))
function getIdx(pp::InferenceVariable, sym::Symbol, i::Int=0)
  return getIdx(pp.dimtype, sym)
end



"""
    $SIGNATURES

Return `::Bool` on whether this variable has been marginalized.
"""
isMarginalized(vert::DFGVariable) = getSolverData(vert).ismargin
isMarginalized(dfg::AbstractDFG, sym::Symbol) = isMarginalized(DFG.getVariable(dfg, sym))

function setThreadModel!(fgl::AbstractDFG;
                         model=IncrementalInference.SingleThreaded )
  #
  for (key, id) in fgl.fIDs
    getSolverData(getFactor(fgl, key)).fnc.threadmodel = model
  end
  nothing
end

"""
    $SIGNATURES

Return bool on whether a certain factor has user defined multihypothesis.

Related

getMultihypoDistribution
"""
isMultihypo(fct::DFGFactor) = isa(getSolverData(fct).fnc.hypotheses, Distribution)

"""
    $SIGNATURES

Return the categorical distributed used for multihypothesis selection in a factor.

Related

isMultihypo
"""
getMultihypoDistribution(fct::DFGFactor) = getSolverData(fct).fnc.hypotheses

"""
    $SIGNATURES

Free all variables from marginalization.
"""
function dontMarginalizeVariablesAll!(fgl::AbstractDFG)
  fgl.solverParams.isfixedlag = false
  fgl.solverParams.qfl = 9999999999
  fgl.solverParams.limitfixeddown = false
  for sym in ls(fgl)
    getSolverData(getVariable(fgl, sym)).ismargin = false
  end
  nothing
end

"""
    $SIGNATURES

Free all variables from marginalization.

Related

dontMarginalizeVariablesAll!
"""
function unfreezeVariablesAll!(fgl::AbstractDFG)
  dontMarginalizeVariablesAll!(fgl)
end

"""
    $SIGNATURES

Reset initialization flag on all variables in `::AbstractDFG`.

Notes
- Numerical values remain, but inference will overwrite since init flags are now `false`.
"""
function resetVariableAllInitializations!(fgl::AbstractDFG)
  vsyms = ls(fgl)
  for sym in vsyms
    setVariableInitialized!(getVariable(fgl, sym), :false)
  end
  nothing
end

"""
    $SIGNATURES

Enable defaults for fixed-lag-like operation by using smart message passing on the tree.

Notes:
- These are only default settings, and can be modified in each use case scenario.
- Default does not update downsolve through to leaves of the tree.
"""
function defaultFixedLagOnTree!(dfg::AbstractDFG,
                                len::Int=30;
                                limitfixeddown::Bool=true )
  #
  getSolverParams(dfg).isfixedlag = true
  getSolverParams(dfg).qfl = len
  getSolverParams(dfg).limitfixeddown = limitfixeddown
  getSolverParams(dfg)
end

"""
    $SIGNATURES

Return `::Tuple` with matching variable ID symbols and `Suggested` PPE values.

Related

getVariablePPE
"""
function getPPESuggestedAll(dfg::AbstractDFG,
                            regexFilter::Union{Nothing, Regex}=nothing )
  #
  # get values
  vsyms = listVariables(dfg, regexFilter) |> sortDFG
  slamPPE = map(x->getVariablePPE(dfg, x).suggested, vsyms)
  # sizes to convert to matrix
  rumax = zeros(Int, 2)
  for ppe in slamPPE
    rumax[2] = length(ppe)
    rumax[1] = maximum(rumax)
  end

  # populate with values
  XYT = zeros(length(slamPPE),rumax[1])
  for i in 1:length(slamPPE)
    XYT[i,1:length(slamPPE[i])] = slamPPE[i]
  end
  return (vsyms, XYT)
end

"""
    $SIGNATURES

Find and return a `::Tuple` of variables and distances to `loc::Vector{<:Real}`.

Related

findVariablesNearTimestamp
"""
function findVariablesNear(dfg::AbstractDFG,
                           loc::Vector{<:Real},
                           regexFilter::Union{Nothing, Regex}=nothing;
                           number::Int=3  )
  #

  xy = getPPESuggestedAll(dfg, regexFilter)
  dist = sum( (xy[2][:,1:length(loc)] .- loc').^2, dims=2) |> vec
  prm = (dist |> sortperm)[1:number]
  return (xy[1][prm], sqrt.(dist[prm]))
end

"""
    $SIGNATURES

Convenience wrapper to `DFG.loadDFG!` taking only one argument, the file name, to load a DFG object in standard format.
"""
loadDFG(filename::AbstractString) = loadDFG!(initfg(), filename)

## ============================================================================
# Starting integration with Manifolds.jl, via ApproxManifoldProducts.jl first
## ============================================================================

# FIXME, much consolidation required here
convert(::Type{<:AMP.Manifold}, ::InstanceType{ContinuousEuclid}) = AMP.Euclid


#
