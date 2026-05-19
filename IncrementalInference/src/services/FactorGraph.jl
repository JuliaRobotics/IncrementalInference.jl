
"""
$SIGNATURES

Initialize an empty in-memory DistributedFactorGraph `::DistributedFactorGraph` object.
"""
function initfg(
  dfg::T = LocalDFG(; solverParams = SolverParams());
  sessionname = "NA",
  robotname = "",
  username = "",
  cloudgraph = nothing,
) where {T <: AbstractDFG}
  #
  #
  return dfg
end

#init an empty fg with a provided type and SolverParams
function initfg(
  ::Type{T};
  solverParams = SolverParams(),
  sessionname = "NA",
  robotname = "",
  username = "",
  cloudgraph = nothing,
) where {T <: AbstractDFG}
  #
  return T(; solverParams = solverParams)
end

function initfg(
  ::Type{T},
  solverParams::S;
  sessionname = "NA",
  robotname = "",
  username = "",
  cloudgraph = nothing,
) where {T <: AbstractDFG, S <: SolverParams}
  #
  return T{S}(; solverParams = solverParams)
end

# Should deprecate in favor of TensorCast.jl
reshapeVec2Mat(vec::Vector, rows::Int) = reshape(vec, rows, round(Int, length(vec) / rows))

## ==============================================================================================
## MOVE TO / CONSOLIDATE WITH DFG
## ==============================================================================================

"""
    $(SIGNATURES)

Fetch the variable marginal joint sampled points.  Use [`getBelief`](@ref) to retrieve the full Belief object.
"""
#FIXME replace with refPoints
getVal(v::VariableCompute; solveKey::Symbol = :default) = DFG.refPoints(getState(v, solveKey))
function getVal(v::VariableCompute, idx::Int; solveKey::Symbol = :default)
  return DFG.refPoints(getState(v, solveKey))[idx]
end
getVal(vnd::State) = DFG.refPoints(vnd)
getVal(vnd::State, idx::Int) = DFG.refPoints(vnd)[idx]
function getVal(dfg::AbstractDFG, lbl::Symbol; solveKey::Symbol = :default)
  return DFG.refPoints(getVariable(dfg, lbl).states[solveKey])
end

"""
    $(SIGNATURES)

Get the number of points used for the current marginal belief estimate represtation for a particular variable in the factor graph.
"""
function getNumPts(v::VariableCompute; solveKey::Symbol = :default)::Int
  return length(getVal(getState(v, solveKey)))
end

function AMP.getBW(vnd::State)
  return DFG.refBandwidth(vnd)
end

# setVal! assumes you will update values to database separate, this used for local graph mods only
function getBWVal(v::VariableCompute; solveKey::Symbol = :default)
  return DFG.refBandwidth(getState(v, solveKey))
end
function setBW!(vd::State, bw::Array{Float64, 2}; solveKey::Symbol = :default)
  DFG.refBandwidth(vd) .= bw
  return nothing
end
function setBW!(v::VariableCompute, bw::Array{Float64, 2}; solveKey::Symbol = :default)
  setBW!(getState(v, solveKey), bw)
  return nothing
end

function setVal!(
  vd::State, 
  val::AbstractVector{P}; 
  observability::AbstractVector{<:Real} = [0.0;]
) where {P}
  points = DFG.refPoints(vd)
  resize!(points, length(val))
  points .= val

  observability = DFG.refObservability(vd)
  resize!(observability, length(observability))
  observability .= observability
  return nothing
end
function setVal!(
  v::VariableCompute,
  val::AbstractVector{P};
  solveKey::Symbol = :default,
  observability::AbstractVector{<:Real} = [0.0;],
) where {P}
  setVal!(getState(v, solveKey), val; observability)
  return nothing
end
function setVal!(
  vd::State,
  val::AbstractVector{P},
  bw::AbstractMatrix{Float64};
  observability::AbstractVector{<:Real} = [0.0;],
) where {P}
  setVal!(vd, val; observability)
  setBW!(vd, bw)
  return nothing
end
function setVal!(
  v::VariableCompute,
  val::AbstractVector{P},
  bw::AbstractMatrix{Float64};
  solveKey::Symbol = :default,
  observability::AbstractVector{<:Real} = [0.0;],
) where {P}
  setVal!(v, val; solveKey, observability)
  setBW!(v, bw; solveKey)
  return nothing
end
function setVal!(
  vd::State,
  val::AbstractVector{P},
  bw::AbstractVector{Float64};
  observability::AbstractVector{<:Real} = [0.0;],
) where {P}
  setVal!(vd, val, reshape(bw, length(bw), 1); observability)
  return nothing
end
function setVal!(
  v::VariableCompute,
  val::AbstractVector{P},
  bw::AbstractVector{Float64};
  solveKey::Symbol = :default,
  observability::AbstractVector{<:Real} = [0.0;],
) where {P}
  setVal!(getState(v, solveKey), val, bw; observability)
  return nothing
end
function setVal!(
  dfg::AbstractDFG,
  sym::Symbol,
  val::AbstractVector{P};
  solveKey::Symbol = :default,
  observability::AbstractVector{<:Real} = [0.0;],
) where {P}
  return setVal!(getVariable(dfg, sym), val; solveKey, observability)
end

"""
    $SIGNATURES

Set the point centers and bandwidth parameters of a variable node, also set `isInitialized=true` if `setinit::Bool=true` (as per default).

Notes
- `initialized` is used for initial solve of factor graph where variables are not yet initialized.
- `inferdim` is used to identify if the initialized was only partial.
"""
function setValKDE!(
  vd::State,
  pts::AbstractVector{P},
  bws::Vector{Float64},
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;],
) where {P}
  #

  setVal!(vd, pts, bws; observability = ipc) # BUG ...al!(., val, . ) ## TODO -- this can be a little faster
  setinit ? (vd.initialized = true) : nothing
  # vd.observability = ipc # TODO, state.belief.observability = ipc instead
  return nothing
end

function setValKDE!(
  vd::State,
  val::AbstractVector{P},
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;],
) where {P}
  # recover variableType information
  varType = getStateKind(vd)
  p = AMP.manikde!(varType, val)
  setValKDE!(vd, p, setinit, ipc)
  return nothing
end

function setValKDE!(
  v::VariableCompute,
  val::AbstractVector{P},
  bws::Array{<:Real, 2},
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;];
  solveKey::Symbol = :default,
) where {P}
  # recover variableType information
  setValKDE!(getState(v, solveKey), val, bws[:, 1], setinit, ipc)

  return nothing
end

function setValKDE!(
  v::VariableCompute,
  val::AbstractVector{P},
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;];
  solveKey::Symbol = :default,
) where {P}
  vnd = getState(v, solveKey)
  # recover variableType information
  setValKDE!(vnd, val, setinit, ipc)
  return nothing
end
function setValKDE!(
  v::VariableCompute,
  em::TreeBelief,
  setinit::Bool = true;
  # inferdim::Union{Float32, Float64, Int32, Int64}=0;
  solveKey::Symbol = :default,
)
  #
  setValKDE!(v, em.val, em.bw, setinit, em.infoPerCoord; solveKey = solveKey)
  return nothing
end
function setValKDE!(
  v::VariableCompute,
  mkd::ManifoldKernelDensity,
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;];
  solveKey::Symbol = :default,
)
  #
  # @error("TESTING setValKDE! ", solveKey, string(listStates(v)))
  setValKDE!(getState(v, solveKey), mkd, setinit, Float64.(ipc))
  return nothing
end
function setValKDE!(
  dfg::AbstractDFG,
  sym::Symbol,
  mkd::ManifoldKernelDensity,
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;];
  solveKey::Symbol = :default,
)
  #
  setValKDE!(getVariable(dfg, sym), mkd, setinit, ipc; solveKey = solveKey)
  return nothing
end

function setValKDE!(
  vnd::State,
  mkd::ManifoldKernelDensity{M, B, Nothing}, # TBD dispatch without partial?
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;],
) where {M, B}
  #
  # L==Nothing means no partials
  ptsArr = AMP.getPoints(mkd) # , false) # for not partial
  # also set the bandwidth
  bws = getBW(mkd)[:, 1]
  setValKDE!(vnd, ptsArr, bws, setinit, ipc)
  return nothing
end

function setValKDE!(
  vnd::State,
  mkd::ManifoldKernelDensity{M, B, L},
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;],
) where {M, B, L <: AbstractVector}
  #
  oldBel = getBelief(vnd)

  # New infomation might be partial
  newBel = replace(oldBel, mkd)

  # Set partial dims as Manifold points
  ptsArr = AMP.getPoints(newBel, false)

  # also get the bandwidth
  bws = getBandwidth(newBel, false)

  # update values in graph
  setValKDE!(vnd, ptsArr, bws, setinit, ipc)
  return nothing
end

function setBelief!(
  vari::VariableCompute, 
  bel::ManifoldKernelDensity, 
  setinit::Bool=true, 
  ipc::AbstractVector{<:Real}=[0.0;];
  solveKey::Symbol = :default
)
  setValKDE!(vari, bel, setinit, ipc; solveKey)
  # setValKDE!(vari,getPoints(bel, false), setinit, ipc)
end

"""
    $SIGNATURES

Set variable initialized status.
"""
function setVariableInitialized!(varid::State, status::Bool)
  return varid.initialized = status
end

function setVariableInitialized!(vari::VariableCompute, solveKey::Symbol, status::Bool)
  return setVariableInitialized!(getState(vari, solveKey), status)
end

"""
    $SIGNATURES

Set method for the inferred dimension value in a variable.
"""
setIPC!(varid::State, val::AbstractVector{<:Real}) = varid.observability = val
function setIPC!(
  vari::VariableCompute,
  val::AbstractVector{<:Real},
  solveKey::Symbol = :default,
)
  return setVariableIPC!(getState(vari, solveKey), val)
end

## ==============================================================================================
## ==============================================================================================

"""
    $(SIGNATURES)

Get a ManifoldKernelDensity estimate from variable node data.
"""
function getBelief(vnd::State)
  return manikde!(getManifold(getStateKind(vnd)), getVal(vnd); bw = getBW(vnd)[:, 1])
end

function getBelief(v::VariableCompute, solvekey::Symbol = :default)
  return getBelief(getState(v, solvekey))
end
function getBelief(dfg::AbstractDFG, lbl::Symbol, solvekey::Symbol = :default)
  return getBelief(getVariable(dfg, lbl), solvekey)
end

"""
    $SIGNATURES

Reset the solve state of a variable to uninitialized/unsolved state.
"""
function resetVariable!(varid::State)
  #
  val = getBelief(varid)
  pts = AMP.getPoints(val)
  # TODO not all manifolds will initialize to zero
  for pt in pts
    fill!(pt, 0.0)
  end
  pn = manikde!(getManifold(varid), pts; bw = zeros(Ndim(val)))
  setValKDE!(varid, pn, false, [0.0;])
  # setVariableInferDim!(varid, 0)
  # setVariableInitialized!(vari, false)
  return nothing
end

function resetVariable!(vari::VariableCompute, solveKey::Symbol = :default)
  return resetVariable!(getState(vari, solveKey))
end

function resetVariable!(dfg::AbstractDFG, sym::Symbol, solveKey::Symbol = :default)
  return resetVariable!(getState(dfg, sym, solveKey))
end

# DefaultNodeDataParametric is deprecated — logic moved to prepareState!(v, NLLSSolver(), statelabel)

"""
    $SIGNATURES

Makes and sets a parametric `State` object (`.solverData`).

!!! warning "Deprecated"
    Use `prepareState!(v, NLLSSolver(), statelabel)` instead.
"""
function setDefaultNodeDataParametric!(
  v::VariableCompute,
  variableType::StateType;
  solveKey::Symbol = :parametric,
  kwargs...,
)
  Base.depwarn("`setDefaultNodeDataParametric!` is deprecated, use `prepareState!(v, NLLSSolver(), solveKey)` instead.", :setDefaultNodeDataParametric!)
  prepareState!(v, NLLSSolver(), solveKey)
  nothing
end

"""
    $SIGNATURES

Create new solverData.

!!! warning "Deprecated"
    Use `prepareState!(v, NPBPSolver(), statelabel; num_kernels=N)` instead.
"""
function setDefaultNodeData!(
  v::VariableCompute,
  dodims::Int,
  N::Int;
  solveKey::Symbol = :default,
  gt = Dict(),
  initialized::Bool = true,
  varType = nothing,
)
  Base.depwarn("`setDefaultNodeData!` is deprecated, use `prepareState!(v, NPBPSolver(), solveKey; num_kernels=N)` instead.", :setDefaultNodeData!)
  prepareState!(v, NPBPSolver(), solveKey; num_kernels=N)
  return nothing
end
# if size(initval,2) < N && size(initval, 1) == dims
#   @warn "setDefaultNodeData! -- deprecated use of stdev."
#   p = manikde!(varType.manifold, initval,diag(stdev));
#   pN = resample(p,N)
# if size(initval,2) < N && size(initval, 1) != dims
# @info "Node value memory allocated but not initialized"
# else
#   pN = manikde!(varType.manifold, initval)
# end
# dims = size(initval,1) # rows indicate dimensions

"""
    $SIGNATURES

Reference data can be stored in the factor graph as a super-solve.

Notes
- Intended as a mechanism to store reference data alongside the numerical computations.
"""
function setVariableRefence!(
  dfg::AbstractDFG,
  sym::Symbol,
  val::AbstractVector;
  refKey::Symbol = :reference,
)
  #
  # which variable to update
  var = getVariable(dfg, sym)

  # Construct an empty VND object
  vnd = State(
    val,
    zeros(getDimension(var), 1),
    Symbol[],
    Int[0;],
    getDimension(var),
    false,
    :_null,
    Symbol[],
    getStateKind(var),
    true,
    zeros(getDimension(var)),
    false,
    true,
  )
  #
  # set the value in the VariableCompute
  return mergeState!(var, vnd)
end

# get instance from variableType
_variableType(varType::StateType) = varType
_variableType(varType::Type{<:StateType}) = varType()

## ==================================================================================================
## DFG Overloads on addVariable! and addFactor!
## ==================================================================================================

#TODO move to DFG after deprecation of IIF kwargs
function DFG.addVariable!(
  dfg::AbstractDFG,
  label::Symbol,
  statekind::Union{T, Type{T}};
  # deprecated in v0.37
  N::Union{Int, Nothing} = nothing,
  initsolvekeys::Union{Vector{Symbol}, Nothing} = nothing,
  #
  kwargs...,  
) where {T <: StateType}
  
  if !isnothing(N)
    Base.depwarn(
      "`addVariable!(dfg, lbl, T; N=$N)` is deprecated. " *
      "Particle count `N` is now a solver option passed to `solveTree!` or `initAll!`.",
      :addVariable!,
    )
    num_kernels = N
  else
    num_kernels = 100 #FIXME 
  end

  if !isnothing(initsolvekeys)
    Base.depwarn(
      "`addVariable!(dfg, lbl, T; initsolvekeys=...)` is deprecated. " *
      "Solver data is created by `initAll!`/`solveTree!`. The kwarg is ignored.",
      :addVariable!,
    )
  end

  v = VariableDFG(
    label,
    statekind;
    kwargs...,
  )

  #FIXME needs better refactoring of spaghetti between addVariable! and prepareState! 
  #FIXME so hard coding :default state creation here for now. Ideal is complete decoupling of solver options from graph operations.
  #also num_kernels set above
  prepareState!(v, NPBPSolver(), :default; num_kernels)

  return addVariable!(dfg, v)
end

function parseusermultihypo!(multihypo::Vector{Float64})
  if isempty(multihypo)
    return nothing
  end
  multihypo2 = multihypo
  multihypo2[1 - 1e-10 .< multihypo] .= 0.0
  # check that terms sum to full probability
  @assert abs(sum(multihypo2) % 1) < 1e-10 || 1 - 1e-10 < sum(multihypo2) % 1 "ensure multihypo sums to a (or nearly, 1e-10) interger, see #1086"
  # check that only one variable broken into fractions
  @assert sum(multihypo2[1e-10 .< multihypo2]) ≈ 1
  # force normalize something that is now known to be close
  multihypo2 ./= sum(multihypo2)
  return Categorical(Float64[multihypo2...])
end

# return a BitVector masking the fractional portion, assuming converted 0's on 100% confident variables 
function _getFractionalVars(varList::Union{<:Tuple, <:AbstractVector}, mh::Nothing)
  return zeros(length(varList)) .== 1
end
_getFractionalVars(varList::Union{<:Tuple, <:AbstractVector}, mh::Categorical) = 0 .< mh.p

function _selectHypoVariables(
  allVars::Union{<:Tuple, <:AbstractVector},
  mh::Categorical,
  sel::Integer = rand(mh),
)
  #
  mask = mh.p .≈ 0.0
  mask[sel] = true
  return (1:length(allVars))[mask]
end

function _selectHypoVariables(
  allVars::Union{<:Tuple, <:AbstractVector},
  mh::Nothing,
  sel::Integer = 0,
)
  return collect(1:length(allVars))
end

"""
    $SIGNATURES

Overload for specific factor preamble usage.

Notes:
- See https://github.com/JuliaRobotics/IncrementalInference.jl/issues/1462

DevNotes
- Integrate into CalcFactor
  - Add threading

Example:

```julia
import IncrementalInference: preambleCache

preambleCache(dfg::AbstractDFG, vars::AbstractVector{<:VariableCompute}, usrfnc::MyFactor) = MyFactorCache(randn(10))

# continue regular use, e.g.
mfc = MyFactor(...)
addFactor!(fg, [:a;:b], mfc)
# ... 
```
"""
function preambleCache(
  dfg::AbstractDFG,
  vars::AbstractVector{<:VariableCompute},
  observation::AbstractObservation,
)
  return nothing
end

# TODO perhaps consolidate with constructor?
"""
$SIGNATURES

Generate the default factor data for a new FactorCompute.
"""
function prepareFactorCache!(
  dfg::AbstractDFG,
  factor::FactorCompute,
  neighbors::AbstractVector{<:VariableCompute} = collect(getVariable.(dfg, getVariableOrder(factor)));
  _blockRecursion::Bool = true,
  attemptGradients::Bool = false,
  keepCalcFactor::Bool = false,
)
  obs = DFG.getObservation(factor)
  multihypo = parseusermultihypo!(factor.hyper.multihypo)
  userCache = preambleCache(dfg, neighbors, obs)   # TODO prepareCache(...)
  ccw =  _createCCW(
    neighbors,
    obs;
    multihypo,
    nullhypo = factor.hyper.nullhypo,
    inflation = factor.hyper.inflation,
    attemptGradients,
    _blockRecursion,
    userCache,
    keepCalcFactor,
  )
  DFG.setCache!(factor, ccw)
  return ccw
end

"""
    $SIGNATURES

Return `::Bool` on whether at least one hypothesis is available for intended computations (assuming direction `sfidx`).
"""
function isLeastOneHypoAvailable(
  sfidx::Int,
  certainidx::Vector{Int},
  uncertnidx::Vector{Int},
  isinit::Vector{Bool},
)
  #
  # @show isinit
  # @show sfidx in certainidx, sum(isinit[uncertnidx])
  # @show sfidx in uncertnidx, sum(isinit[certainidx])
  return sfidx in certainidx && 0 < sum(isinit[uncertnidx]) ||
         sfidx in uncertnidx && sum(isinit[certainidx]) == length(certainidx)
end

function assembleFactorName(dfg::AbstractDFG, Xi::Vector{<:VariableCompute}; maxincidence::Int = 500)
  #

  existingFactorLabels = listFactors(dfg)
  existingFactorLabelDict = Dict(existingFactorLabels .=> existingFactorLabels)
  namestring = ""
  for vert in Xi #f.Xi
    namestring = string(namestring, vert.label)
  end
  for i = 1:maxincidence
    tempnm = string(namestring, "f$i")
    if !haskey(existingFactorLabelDict, Symbol(tempnm))
      namestring = tempnm
      break
    end
    if i != maxincidence
      nothing
    else
      error(
      "Artificial restriction to not connect more than $(maxincidence) factors to a variable (bad for sparsity).",
    )
    end
  end
  return Symbol(namestring)
end

"""
    $(SIGNATURES)

Add factor with user defined type `<:AbstractObservation` to the factor graph object.

Notes
- This is a graph operation. No automatic variable initialization is performed.
- Use `initAll!` or `solveTree!` to initialize variables before solving.
- `graphinit` kwarg is deprecated and ignored. Call `doautoinit!` explicitly if needed.

Experimental
- `inflation`, to better disperse kernels before convolution solve, see IIF #1051.
"""
function DFG.addFactor!(
  dfg::AbstractDFG,
  Xi::AbstractVector{<:VariableCompute},
  observation::AbstractObservation;
  namestring::Symbol = assembleFactorName(dfg, Xi),
  #TODO solver parameters/options follows
  inflation::Real = 5.0,
  _blockRecursion::Bool = true,
  keepCalcFactor::Bool = false,
  # Deprecated in v0.37
  graphinit::Union{Bool, Nothing} = nothing,
  kwargs...
)
  #

  if !isnothing(graphinit)
    Base.depwarn(
      "`addFactor!(dfg, vars, fct; graphinit=$graphinit)` is deprecated. " *
      "Variable initialization is no longer done in addFactor!. " *
      "Use `initAll!(dfg)` or `doautoinit!(dfg, vars)` explicitly.",
      :addFactor!,
    )
  end

  variableorder = Symbol[v.label for v in Xi]
  #
  newFactor = FactorDFG(
    variableorder,
    observation;
    label = namestring,
    inflation,
    kwargs...
  )

  return addFactor!(dfg, newFactor)
end

function DFG.addFactor!(
  dfg::AbstractDFG,
  vlbs::AbstractVector{Symbol},
  usrfnc::AbstractObservation;
  kw...,
)
  variables = map(vid -> getVariable(dfg, vid), vlbs)
  return addFactor!(dfg, variables, usrfnc; kw...)
end

#
