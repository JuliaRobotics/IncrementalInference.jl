# Factor Graph OS type utilities
#  IIF methods should direclty detect extended types from user import
# of convert in their namespace


# FIXME, upgrade to AMP instead
KDE.getPoints(dfg::AbstractDFG, lbl::Symbol) = getBelief(dfg, lbl) |> getPoints

clampStringLength(st::AbstractString, len::Int=5) = st[1:minimum([len; length(st)])]

function clampBufferString(st::AbstractString, max::Int, len::Int=minimum([max,length(st)]))
  @assert 0 <= max "max must be greater or equal to zero"
  st = clampStringLength(st, len)
  for i in len:max-1  st *= " "; end
  return st
end

"""
    $SIGNATURES

Extract contiguous string of numbers at end of a label`::Symbol` -- e.g. `:x45_4` --> "4".  
Returns `(string, suffix_substring)`

Related

[`incrSuffix`](@ref)
"""
function _getSuffix(lbl::Symbol; pattern::Regex=r"\d+")
  slbl = string(lbl)
  phrase_ = slbl |> reverse |> x->match(pattern,x).match 
  slbl, reverse(phrase_)
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
function incrSuffix(lbl::Symbol, val::Integer=+1; pattern::Regex=r"\d+")
  slbl, phrase = _getSuffix(lbl, pattern=pattern)
  nint = phrase |> x->(parse(Int,x)+val)
  prefix = slbl[1:(end-length(phrase))]
  Symbol(prefix, nint)
end

"""
    $SIGNATURES
Get the CommonConvWrapper for this factor.
"""
_getCCW(gfnd::GenericFunctionNodeData) = gfnd.fnc
_getCCW(fct::DFGFactor) = getSolverData(fct) |> _getCCW
_getCCW(dfg::AbstractDFG, lbl::Symbol) = getFactor(dfg, lbl) |> _getCCW


DFG.getFactorType(ccw::CommonConvWrapper) = ccw.usrfnc!

_getZDim(ccw::CommonConvWrapper) = ccw.zDim
# TODO is MsgPrior piggy backing zdim on inferdim???
_getZDim(ccw::CommonConvWrapper{<:MsgPrior}) = length(ccw.usrfnc!.infoPerCoord) # ccw.usrfnc!.inferdim

_getZDim(fcd::GenericFunctionNodeData) = _getCCW(fcd) |> _getZDim
_getZDim(fct::DFGFactor) = _getCCW(fct) |> _getZDim

DFG.getDimension(fct::GenericFunctionNodeData) = _getZDim(fct)
DFG.getDimension(fct::DFGFactor) = _getZDim(fct)


"""
    $SIGNATURES

Return the manifold on which this ManifoldKernelDensity is defined.

DevNotes
- TODO currently ignores the .partial aspect (captured in parameter `L`)
"""
getManifold(mkd::ManifoldKernelDensity{M,B,Nothing}, asPartial::Bool=false) where {M,B} = mkd.manifold
getManifold(mkd::ManifoldKernelDensity{M,B,L}, asPartial::Bool=false) where {M,B,L <: AbstractVector} = asPartial ? mkd.manifold : getManifoldPartial(mkd.manifold, mkd._partial)


"""
    $TYPEDSIGNATURES

Return the number of dimensions this factor vertex `fc` influences.
"""
getFactorDim(w...) = getDimension(w...)
# getFactorDim(fcd::GenericFunctionNodeData) = isa(_getCCW(fcd).usrfnc!, MsgPrior) ? _getCCW(fcd).usrfnc!.inferdim : Int(_getCCW(fcd).zDim)
# getFactorDim(fc::DFGFactor) = getFactorDim(getSolverData(fc))
getFactorDim(fg::AbstractDFG, fctid::Symbol) = getFactorDim(getFactor(fg, fctid))



# function _getDimensionsPartial(ccw::CommonConvWrapper)
#   # @warn "_getDimensionsPartial not ready for use yet"
#   ccw.partialDims
# end
# _getDimensionsPartial(data::GenericFunctionNodeData) = _getCCW(data) |> _getDimensionsPartial
# _getDimensionsPartial(fct::DFGFactor) = _getDimensionsPartial(_getCCW(fct))
# _getDimensionsPartial(fg::AbstractDFG, lbl::Symbol) = _getDimensionsPartial(getFactor(fg, lbl))


# """
#     $SIGNATURES
# Get `.factormetadata` for each CPT in CCW for a specific factor in `fg`. 
# """
# _getFMdThread(ccw::CommonConvWrapper, 
#               thrid::Int=Threads.threadid()) = ccw.cpt[thrid].factormetadata
# #
# _getFMdThread(fc::Union{GenericFunctionNodeData,DFGFactor}, 
#               thrid::Int=Threads.threadid()) = _getFMdThread(_getCCW(fc), thrid)
# #
# _getFMdThread(dfg::AbstractDFG,
#               lbl::Symbol,
#               thrid::Int=Threads.threadid()) = _getFMdThread(_getCCW(dfg, lbl), thrid)
# #



# extend convenience function (Matrix or Vector{P})
function manikde!(variableType::Union{InstanceType{<:InferenceVariable}, InstanceType{<:AbstractFactor}},
                  pts::AbstractVector{P};
                  kw... ) where {P <: Union{<:AbstractArray,<:Number,<:ProductRepr,<:Manifolds.ArrayPartition} }
  #
  M = getManifold(variableType)
  infoPerCoord=ones(AMP.getNumberCoords(M, pts[1]))
  return AMP.manikde!(M, pts; infoPerCoord, kw...)
end


"""
    $SIGNATURES

Return params.N measurement samples for a factor in `<:AbstractDFG`.
"""
getMeasurements(dfg::AbstractDFG, fsym::Symbol, N::Int=getSolverParams(dfg).N) = sampleFactor(dfg, fsym, N)


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
  
  # check that the variable to freeze exists fix issue #966
  filter!(v->exists(dfg,v), tofreeze)

  if length(tofreeze) == 0
      @info "[fifoFreeze] QFL - no nodes to freeze."
      return nothing
  end
  @info "[fifoFreeze] QFL - Freezing nodes $(tofreeze[1]) -> $(tofreeze[end])."
  setfreeze!(dfg, tofreeze)
  nothing
end


DFG.getPoint(typ::InferenceVariable, w...;kw...) = getPoint(typeof(typ), w...;kw...)
DFG.getCoordinates(typ::InferenceVariable, w...;kw...) = getCoordinates(typeof(typ), w...;kw...)

# WIP
# _getMeasurementRepresentation(::AbstractPrior, coord::AbstractVector{<:Number}) = 


"""
    $SIGNATURES

Get the ParametricPointEstimates---based on full marginal belief estimates---of a variable in the distributed factor graph.

DevNotes
- TODO update for manifold subgroups.
- TODO standardize after AMP3D

Related

[`getVariablePPE`](@ref), [`setVariablePosteriorEstimates!`](@ref), [`getVariablePPE!`](@ref)
"""
function calcPPE( var::DFGVariable,
                  varType::InferenceVariable=getVariableType(var);
                  ppeType::Type{<:MeanMaxPPE}=MeanMaxPPE,
                  solveKey::Symbol=:default,
                  ppeKey::Symbol=solveKey,
                  timestamp=now()  )
  #
  P = getBelief(var, solveKey)
  maniDef = convert(MB.AbstractManifold, varType)
  manis = convert(Tuple, maniDef) # LEGACY, TODO REMOVE
  ops = buildHybridManifoldCallbacks(manis)
  Pme = calcMean(P)  # getKDEMean(P) #, addop=ops[1], diffop=ops[2]

  # returns coordinates at identify
  Pma = getKDEMax(P, addop=ops[1], diffop=ops[2])
  # calculate point

  ## TODO make PPE only use getCoordinates for now (IIF v0.25)
  Pme_ = getCoordinates(varType,Pme)
  # Pma_ = getCoordinates(M,Pme)
  
  # suggested, max, mean, current time
  # TODO, poor constructor argument assumptions on `ppeType`
  ppeType(ppeKey, Pme_, Pma, Pme_, timestamp)
end


# calcPPE(var::DFGVariable; method::Type{<:AbstractPointParametricEst}=MeanMaxPPE, solveKey::Symbol=:default) = calcPPE(var, getVariableType(var), method=method, solveKey=solveKey)

"""
    $TYPEDSIGNATURES

Calculate new Parametric Point Estimates for a given variable.

Notes
- Different methods are possible, currently [`MeanMaxPPE`](@ref) `<: AbstractPointParametricEst`.

Aliases
- `calcVariablePPE`

Related

[`setPPE!`](@ref)
"""
function calcPPE( dfg::AbstractDFG,
                  label::Symbol;
                  solveKey::Symbol=:default,
                  ppeType::Type{<:AbstractPointParametricEst}=MeanMaxPPE )
  #
  var = getVariable(dfg, label)
  calcPPE(var, getVariableType(var), ppeType=ppeType, solveKey=solveKey)
end

const calcVariablePPE = calcPPE


function setThreadModel!( fgl::AbstractDFG;
                          model=IIF.SingleThreaded )
  #
  for (key, id) in fgl.fIDs
    _getCCW(fgl, key).threadmodel = model
  end
  nothing
end

"""
    $SIGNATURES

Return bool on whether a certain factor has user defined multihypothesis.

Related

getMultihypoDistribution
"""
isMultihypo(fct::DFGFactor) = isa(_getCCW(fct).hypotheses, Distribution)

"""
    $SIGNATURES

Return the categorical distributed used for multihypothesis selection in a factor.

Related

isMultihypo
"""
getMultihypoDistribution(fct::DFGFactor) = _getCCW(fct).hypotheses




"""
    $SIGNATURES

Free all variables from marginalization.
"""
function dontMarginalizeVariablesAll!(fgl::AbstractDFG)
  fgl.solverParams.isfixedlag = false
  fgl.solverParams.qfl = (2^(Sys.WORD_SIZE-1)-1)
  fgl.solverParams.limitfixeddown = false
  for sym in ls(fgl)
    setMarginalized!(fgl, sym, false)
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
function findVariablesNear( dfg::AbstractDFG,
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



"""
    $SIGNATURES

Find all factors that go `from` variable to any other complete variable set within `between`.

Notes
- Developed for downsolve in CSM, expanding the cliqSubFg to include all frontal factors.
"""
function findFactorsBetweenFrom(dfg::G,
                                between::Vector{Symbol},
                                from::Symbol ) where {G <: AbstractDFG}
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
function getFactorsAmongVariablesOnly(dfg::G,
                                      varlist::Vector{Symbol};
                                      unused::Bool=true  ) where G <: AbstractDFG
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
      if !getSolverData(vert).potentialused
        push!(almostfcts, fct)
      end
    end
  else
    almostfcts = prefcts
  end

  # Select factors that have all variables in this clique var list
  usefcts = Symbol[]
  for fct in almostfcts
    if length(setdiff(DFG.getNeighbors(dfg, fct), varlist)) == 0
      push!(usefcts, fct)
    end
  end

  return usefcts
end


"""
    $SIGNATURES

Calculate new and then set PPE estimates for variable from some distributed factor graph.

DevNotes
- TODO solve key might be needed if one only wants to update one
- TODO consider a more fiting name.
- guess it would make sense that :default=>variableNodeData, goes with :default=>MeanMaxPPE

Aliases
- `setVariablePosteriorEstimates!`

DevNotes:

JT - TODO if subfg is in the cloud or from another fg it has to be updated
it feels like a waste to update the whole variable for one field.
currently i could find mergeUpdateVariableSolverData()
might be handy to use a setter such as updatePointParametricEst(dfg, variable, solverkey)
This might also not be the correct place, if it is uncomment:
````
if (subfg <: InMemoryDFGTypes)
  updateVariable!(subfg, var)
end
```

Related

[`calcPPE`](@ref), getVariablePPE, (updatePPE! ?)
"""
function setPPE!( variable::DFGVariable,
                  solveKey::Symbol = :default,
                  ppeType::Type{T} = MeanMaxPPE,
                  newPPEVal::T = calcPPE(variable, ppeType=ppeType, solveKey=solveKey) ) where {T <: AbstractPointParametricEst}
  #
  # vnd = getSolverData(variable, solveKey)

  #TODO in the future one can perhaps populate other solver data types here by looking at the typeof ppeDict entries
  getPPEDict(variable)[solveKey] = newPPEVal

  return variable
end

function setPPE!( subfg::AbstractDFG,
                  label::Symbol,
                  solveKey::Symbol = :default,
                  ppeType::Type{T} = MeanMaxPPE,
                  newPPEVal::NothingUnion{T} = nothing )  where {T <: AbstractPointParametricEst}
  #
  variable = getVariable(subfg,label)
  # slight optimization to avoid double variable lookup (should be optimized out during code lowering)
  newppe = newPPEVal !== nothing ? newPPEVal : calcPPE(variable, solveKey=solveKey, ppeType=ppeType)  
  setPPE!(variable, solveKey, ppeType, newppe)
end


const setVariablePosteriorEstimates! = setPPE!


## ============================================================================
# Starting integration with Manifolds.jl, via ApproxManifoldProducts.jl first
## ============================================================================



"""
    $SIGNATURES
Fetch and unpack JSON dictionary stored as a data blob.
"""
function fetchDataJSON(dfg::AbstractDFG, varsym::Symbol, lbl::Symbol)
  gde,rawData = getData(dfg, varsym, lbl)
  if gde.mimeType == "application/json/octet-stream"
    JSON2.read(IOBuffer(rawData))
  else
    error("Unknown JSON Blob format $(gde.mimeType)")
  end
end



#
