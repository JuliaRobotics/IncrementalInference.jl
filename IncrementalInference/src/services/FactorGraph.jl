
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
    $SIGNATURES

Reset the solve state of a variable to uninitialized/unsolved state.
"""
function resetVariable!(state::State)
  #
  pts = getPoints(state)
  statekind = getStateKind(state)
  # FIXME, wont work for non-group kinds
  ϵ = getPointIdentity(statekind)
  pts_ = [ϵ for i in 1:length(pts)]
  pn = manikde!(statekind, pts_; bw = zeros(Ndim(ϵ)), newbw = false)
  setBelief!(state, pn)
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
