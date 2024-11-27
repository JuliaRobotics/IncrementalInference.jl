
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
getVal(v::DFGVariable; solveKey::Symbol = :default) = v.solverDataDict[solveKey].val
function getVal(v::DFGVariable, idx::Int; solveKey::Symbol = :default)
  return v.solverDataDict[solveKey].val[:, idx]
end
getVal(vnd::VariableNodeData) = vnd.val
getVal(vnd::VariableNodeData, idx::Int) = vnd.val[:, idx]
function getVal(dfg::AbstractDFG, lbl::Symbol; solveKey::Symbol = :default)
  return getVariable(dfg, lbl).solverDataDict[solveKey].val
end

"""
    $(SIGNATURES)

Get the number of points used for the current marginal belief estimate represtation for a particular variable in the factor graph.
"""
function getNumPts(v::DFGVariable; solveKey::Symbol = :default)::Int
  return length(getVal(getSolverData(v, solveKey)))
end

function AMP.getBW(vnd::VariableNodeData)
  return vnd.bw
end

# setVal! assumes you will update values to database separate, this used for local graph mods only
function getBWVal(v::DFGVariable; solveKey::Symbol = :default)
  return getSolverData(v, solveKey).bw
end
function setBW!(vd::VariableNodeData, bw::Array{Float64, 2}; solveKey::Symbol = :default)
  vd.bw = bw
  return nothing
end
function setBW!(v::DFGVariable, bw::Array{Float64, 2}; solveKey::Symbol = :default)
  setBW!(getSolverData(v, solveKey), bw)
  return nothing
end

function setVal!(vd::VariableNodeData, val::AbstractVector{P}) where {P}
  vd.val = val
  return nothing
end
function setVal!(
  v::DFGVariable,
  val::AbstractVector{P};
  solveKey::Symbol = :default,
) where {P}
  setVal!(getSolverData(v, solveKey), val)
  return nothing
end
function setVal!(
  vd::VariableNodeData,
  val::AbstractVector{P},
  bw::AbstractMatrix{Float64},
) where {P}
  setVal!(vd, val)
  setBW!(vd, bw)
  return nothing
end
function setVal!(
  v::DFGVariable,
  val::AbstractVector{P},
  bw::AbstractMatrix{Float64};
  solveKey::Symbol = :default,
) where {P}
  setVal!(v, val; solveKey = solveKey)
  setBW!(v, bw; solveKey = solveKey)
  return nothing
end
function setVal!(
  vd::VariableNodeData,
  val::AbstractVector{P},
  bw::AbstractVector{Float64},
) where {P}
  setVal!(vd, val, reshape(bw, length(bw), 1))
  return nothing
end
function setVal!(
  v::DFGVariable,
  val::AbstractVector{P},
  bw::AbstractVector{Float64};
  solveKey::Symbol = :default,
) where {P}
  setVal!(getSolverData(v, solveKey), val, bw)
  return nothing
end
function setVal!(
  dfg::AbstractDFG,
  sym::Symbol,
  val::AbstractVector{P};
  solveKey::Symbol = :default,
) where {P}
  return setVal!(getVariable(dfg, sym), val; solveKey = solveKey)
end

"""
    $SIGNATURES

Set the point centers and bandwidth parameters of a variable node, also set `isInitialized=true` if `setinit::Bool=true` (as per default).

Notes
- `initialized` is used for initial solve of factor graph where variables are not yet initialized.
- `inferdim` is used to identify if the initialized was only partial.
"""
function setValKDE!(
  vd::VariableNodeData,
  pts::AbstractVector{P},
  bws::Vector{Float64},
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;],
) where {P}
  #

  setVal!(vd, pts, bws) # BUG ...al!(., val, . ) ## TODO -- this can be a little faster
  setinit ? (vd.initialized = true) : nothing
  vd.infoPerCoord = ipc
  return nothing
end

function setValKDE!(
  vd::VariableNodeData,
  val::AbstractVector{P},
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;],
) where {P}
  # recover variableType information
  varType = getVariableType(vd)
  p = AMP.manikde!(varType, val)
  setValKDE!(vd, p, setinit, ipc)
  return nothing
end

function setValKDE!(
  v::DFGVariable,
  val::AbstractVector{P},
  bws::Array{<:Real, 2},
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;];
  solveKey::Symbol = :default,
) where {P}
  # recover variableType information
  setValKDE!(getSolverData(v, solveKey), val, bws[:, 1], setinit, ipc)

  return nothing
end

function setValKDE!(
  v::DFGVariable,
  val::AbstractVector{P},
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;];
  solveKey::Symbol = :default,
  ppeType::Type{T} = MeanMaxPPE,
) where {P, T}
  vnd = getSolverData(v, solveKey)
  # recover variableType information
  setValKDE!(vnd, val, setinit, ipc)
  setPPE!(v; solveKey, ppeType)
  return nothing
end
function setValKDE!(
  v::DFGVariable,
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
  v::DFGVariable,
  mkd::ManifoldKernelDensity,
  setinit::Bool = true,
  ipc::AbstractVector{<:Real} = [0.0;];
  solveKey::Symbol = :default,
)
  #
  # @error("TESTING setValKDE! ", solveKey, string(listSolveKeys(v)))
  setValKDE!(getSolverData(v, solveKey), mkd, setinit, Float64.(ipc))
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
  vnd::VariableNodeData,
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
  vnd::VariableNodeData,
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
  vari::DFGVariable, 
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
function setVariableInitialized!(varid::VariableNodeData, status::Bool)
  #
  return varid.initialized = status
end

function setVariableInitialized!(vari::DFGVariable, status::Bool)
  return setVariableInitialized!(getSolverData(vari), status)
end

"""
    $SIGNATURES

Set method for the inferred dimension value in a variable.
"""
setIPC!(varid::VariableNodeData, val::AbstractVector{<:Real}) = varid.infoPerCoord = val
function setIPC!(
  vari::DFGVariable,
  val::AbstractVector{<:Real},
  solveKey::Symbol = :default,
)
  return setVariableIPC!(getSolverData(vari, solveKey), val)
end

## ==============================================================================================
## ==============================================================================================

"""
    $(SIGNATURES)

Get a ManifoldKernelDensity estimate from variable node data.
"""
function getBelief(vnd::VariableNodeData)
  return manikde!(getManifold(getVariableType(vnd)), getVal(vnd); bw = getBW(vnd)[:, 1])
end

function getBelief(v::DFGVariable, solvekey::Symbol = :default)
  return getBelief(getSolverData(v, solvekey))
end
function getBelief(dfg::AbstractDFG, lbl::Symbol, solvekey::Symbol = :default)
  return getBelief(getVariable(dfg, lbl), solvekey)
end

"""
    $SIGNATURES

Reset the solve state of a variable to uninitialized/unsolved state.
"""
function resetVariable!(varid::VariableNodeData; solveKey::Symbol = :default)::Nothing
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

function resetVariable!(vari::DFGVariable; solveKey::Symbol = :default)
  return resetVariable!(getSolverData(vari); solveKey = solveKey)
end

function resetVariable!(
  dfg::G,
  sym::Symbol;
  solveKey::Symbol = :default,
)::Nothing where {G <: AbstractDFG}
  #
  return resetVariable!(getVariable(dfg, sym); solveKey = solveKey)
end

# return VariableNodeData
function DefaultNodeDataParametric(
  dodims::Int,
  dims::Int,
  variableType::InferenceVariable;
  initialized::Bool = true,
  dontmargin::Bool = false,
  solveKey::Symbol = :parametric
)
  # this should be the only function allocating memory for the node points
  if false && initialized
    error("not implemented yet")
    # pN = AMP.manikde!(variableType.manifold, randn(dims, N));
    #
    # sp = Int[0;] #round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
    # gbw = getBW(pN)[:,1]
    # gbw2 = Array{Float64}(undef, length(gbw),1)
    # gbw2[:,1] = gbw[:]
    # pNpts = getPoints(pN)
    # #initval, stdev
    # return VariableNodeData(pNpts,
    #                         gbw2, Symbol[], sp,
    #                         dims, false, :_null, Symbol[], variableType, true, 0.0, false, dontmargin)
  else
    # dimIDs = round.(Int, range(dodims; stop = dodims + dims - 1, length = dims))
    ϵ = getPointIdentity(variableType)
    return VariableNodeData(variableType;
      id=nothing,
      val=[ϵ],
      bw=zeros(dims, dims),
      # Symbol[],
      # dimIDs,
      dims,
      # false,
      # :_null,
      # Symbol[],
      initialized=false,
      infoPerCoord=zeros(dims),
      ismargin=false,
      dontmargin,
      # 0,
      # 0,
      solveKey,
    )
  end
end

"""
    $SIGNATURES

Makes and sets a parametric `VariableNodeData` object (`.solverData`).

DevNotes
- TODO assumes parametric solves will always just be under the `solveKey=:parametric`, should be generalized.
"""
function setDefaultNodeDataParametric!(
  v::DFGVariable,
  variableType::InferenceVariable;
  solveKey::Symbol = :parametric,
  kwargs...,
)
  vnd = DefaultNodeDataParametric(0, variableType |> getDimension, variableType; solveKey, kwargs...)
  setSolverData!(v, vnd, solveKey)
  nothing
end

"""
    $SIGNATURES

Create new solverData.

Notes
- Used during creation of new variable, as well as in CSM unique `solveKey`.
"""
function setDefaultNodeData!(
  v::DFGVariable,
  dodims::Int,
  N::Int,
  dims::Int=getDimension(v);
  solveKey::Symbol = :default,
  gt = Dict(),
  initialized::Bool = true,
  dontmargin::Bool = false,
  varType = nothing,
)
  #
  # TODO review and refactor this function, exists as legacy from pre-v0.3.0
  # this should be the only function allocating memory for the node points (unless number of points are changed)
  data = nothing
  isinit = false
  sp = Int[0;]
  (val, bw) = if initialized
    pN = resample(getBelief(v))
    bw = getBW(pN)[:, 1:1]
    pNpts = getPoints(pN)
    isinit = true
    (pNpts, bw)
  else
    sp = round.(Int, range(dodims; stop = dodims + dims - 1, length = dims))
    @assert getPointType(varType) != DataType "cannot add manifold point type $(getPointType(varType)), make sure the identity element argument in @defVariable $varType arguments is correct"
    val = Vector{getPointType(varType)}(undef, N)
    for i = 1:length(val)
      val[i] = getPointIdentity(varType)
    end
    bw = zeros(dims, 1)
    #
    (val, bw)
  end
  # make and set the new solverData
  setSolverData!(
    v,
    VariableNodeData(varType;
      id=nothing,
      val,
      bw,
      # Symbol[],
      # sp,
      dims,
      # false,
      # :_null,
      # Symbol[],
      initialized=isinit,
      infoPerCoord=zeros(getDimension(v)),
      ismargin=false,
      dontmargin,
      # 0,
      # 0,
      solveKey,
    ),
    solveKey,
  )
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
  vnd = VariableNodeData(
    val,
    zeros(getDimension(var), 1),
    Symbol[],
    Int[0;],
    getDimension(var),
    false,
    :_null,
    Symbol[],
    getVariableType(var),
    true,
    zeros(getDimension(var)),
    false,
    true,
  )
  #
  # set the value in the DFGVariable
  return setSolverData!(var, vnd, refKey)
end

# get instance from variableType
_variableType(varType::InferenceVariable) = varType
_variableType(varType::Type{<:InferenceVariable}) = varType()

## ==================================================================================================
## DFG Overloads on addVariable! and addFactor!
## ==================================================================================================

"""
$(SIGNATURES)

Add a variable node `label::Symbol` to `dfg::AbstractDFG`, as `varType<:InferenceVariable`.

Notes
-----
- keyword `nanosecondtime` is experimental and intended as the whole subsection portion -- i.e. accurateTime = (timestamp MOD second) + Nanosecond

Example
-------

```julia
fg = initfg()
addVariable!(fg, :x0, Pose2)
```
"""
function addVariable!(
  dfg::AbstractDFG,
  label::Symbol,
  varTypeU::Union{T, Type{T}};
  N::Int = getSolverParams(dfg).N,
  solvable::Int = 1,
  timestamp::Union{DateTime, ZonedDateTime} = now(localzone()),
  nanosecondtime::Union{Nanosecond, Int64, Nothing} = Nanosecond(0),
  dontmargin::Bool = false,
  tags::Vector{Symbol} = Symbol[],
  smalldata = Dict{Symbol, DFG.SmallDataTypes}(),
  checkduplicates::Bool = true,
  initsolvekeys::Vector{Symbol} = getSolverParams(dfg).algorithms,
) where {T <: InferenceVariable}
  #
  varType = _variableType(varTypeU)

  _zonedtime(s::DateTime) = ZonedDateTime(s, localzone())
  _zonedtime(s::ZonedDateTime) = s 

  union!(tags, [:VARIABLE])
  v = DFGVariable(
    label,
    varType;
    tags = Set(tags),
    smallData = smalldata,
    solvable = solvable,
    timestamp = _zonedtime(timestamp),
    nstime = Nanosecond(nanosecondtime),
  )

  (:default in initsolvekeys) && setDefaultNodeData!(
    v,
    0,
    N,
    getDimension(varType);
    initialized = false,
    varType = varType,
    dontmargin = dontmargin,
  ) # dodims

  (:parametric in initsolvekeys) &&
    setDefaultNodeDataParametric!(v, varType; initialized = false, dontmargin = dontmargin)

  return DFG.addVariable!(dfg, v)
end

function parseusermultihypo(multihypo::Nothing, nullhypo::Float64)
  verts = Symbol[]
  mh = nothing
  return mh, nullhypo
end
function parseusermultihypo(multihypo::Vector{Float64}, nullhypo::Float64)
  mh = nothing
  if 0 < length(multihypo)
    multihypo2 = multihypo
    multihypo2[1 - 1e-10 .< multihypo] .= 0.0
    # check that terms sum to full probability
    @assert abs(sum(multihypo2) % 1) < 1e-10 || 1 - 1e-10 < sum(multihypo2) % 1 "ensure multihypo sums to a (or nearly, 1e-10) interger, see #1086"
    # check that only one variable broken into fractions
    @assert sum(multihypo2[1e-10 .< multihypo2]) ≈ 1
    # force normalize something that is now known to be close
    multihypo2 ./= sum(multihypo2)

    mh = Categorical(Float64[multihypo2...])
  end
  return mh, nullhypo
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
import IncrementalInference: preableCache

preableCache(dfg::AbstractDFG, vars::AbstractVector{<:DFGVariable}, usrfnc::MyFactor) = MyFactorCache(randn(10))

# continue regular use, e.g.
mfc = MyFactor(...)
addFactor!(fg, [:a;:b], mfc)
# ... 
```
"""
function preambleCache(
  dfg::AbstractDFG,
  vars::AbstractVector{<:DFGVariable},
  usrfnc::AbstractFactor,
)
  return nothing
end

# TODO perhaps consolidate with constructor?
"""
$SIGNATURES

Generate the default factor data for a new DFGFactor.
"""
function getDefaultFactorData(
  dfg::AbstractDFG,
  Xi::Vector{<:DFGVariable},
  usrfnc::T;
  multihypo::Vector{<:Real} = Float64[],
  nullhypo::Float64 = 0.0,
  # threadmodel = SingleThreaded,
  eliminated::Bool = false,
  potentialused::Bool = false,
  edgeIDs = Int[],
  solveInProgress = 0,
  inflation::Real = getSolverParams(dfg).inflation,
  _blockRecursion::Bool = false,
) where {T <: AbstractFactor}
  #
  
  # prepare multihypo particulars
  # storeMH::Vector{Float64} = multihypo == nothing ? Float64[] : [multihypo...]
  mhcat, nh = parseusermultihypo(multihypo, nullhypo)
  
  # allocate temporary state for convolutional operations (not stored)
  userCache = preambleCache(dfg, Xi, usrfnc)
  ccwl = _createCCW(
    Xi,
    usrfnc;
    multihypo = mhcat,
    nullhypo = nh,
    inflation,
    attemptGradients = getSolverParams(dfg).attemptGradients,
    _blockRecursion,
    userCache,
  )
  
  # and the factor data itself
  return FunctionNodeData{typeof(ccwl)}(
    eliminated,
    potentialused,
    edgeIDs,
    ccwl,
    multihypo,
    ccwl.hyporecipe.certainhypo,
    nullhypo,
    solveInProgress,
    inflation,
  )
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

function assembleFactorName(dfg::AbstractDFG, Xi::Vector{<:DFGVariable})
  #

  existingFactorLabels = listFactors(dfg)
  existingFactorLabelDict = Dict(existingFactorLabels .=> existingFactorLabels)
  namestring = ""
  for vert in Xi #f.Xi
    namestring = string(namestring, vert.label)
  end
  opt = getSolverParams(dfg)
  for i = 1:(opt.maxincidence)
    tempnm = string(namestring, "f$i")
    if !haskey(existingFactorLabelDict, Symbol(tempnm))
      namestring = tempnm
      break
    end
    if i != opt.maxincidence
      nothing
    else
      error(
      "Artificial restriction to not connect more than $(opt.maxincidence) factors to a variable (bad for sparsity), try setting getSolverParams(fg).maxincidence=1000 to adjust this restriction.",
    )
    end
  end
  return Symbol(namestring)
end

"""
    $(SIGNATURES)

Add factor with user defined type `<:AbstractFactor`` to the factor graph
object. Define whether the automatic initialization of variables should be
performed.  Use order sensitive `multihypo` keyword argument to define if any
variables are related to data association uncertainty.

Experimental
- `inflation`, to better disperse kernels before convolution solve, see IIF #1051.
"""
function DFG.addFactor!(
  dfg::AbstractDFG,
  Xi::AbstractVector{<:DFGVariable},
  usrfnc::AbstractFactor;
  multihypo::Vector{Float64} = Float64[],
  nullhypo::Float64 = 0.0,
  solvable::Int = 1,
  tags::Vector{Symbol} = Symbol[],
  timestamp::Union{DateTime, ZonedDateTime} = now(localzone()),
  graphinit::Bool = getSolverParams(dfg).graphinit,
  # threadmodel = SingleThreaded,
  suppressChecks::Bool = false,
  inflation::Real = getSolverParams(dfg).inflation,
  namestring::Symbol = assembleFactorName(dfg, Xi),
  _blockRecursion::Bool = !getSolverParams(dfg).attemptGradients,
)
  #

  @assert (suppressChecks || length(multihypo) === 0 || length(multihypo) == length(Xi)) "When using multihypo=[...], the number of variables and multihypo probabilies must match.  See documentation on how to include fractional data-association uncertainty."

  _zonedtime(s::ZonedDateTime) = s
  _zonedtime(s::DateTime) = ZonedDateTime(s, localzone())

  varOrderLabels = Symbol[v.label for v in Xi]
  solverData = getDefaultFactorData(
    dfg,
    Xi,
    deepcopy(usrfnc);
    multihypo,
    nullhypo,
    # threadmodel,
    inflation,
    _blockRecursion,
  )
  #
  newFactor = DFGFactor(
    Symbol(namestring),
    varOrderLabels,
    solverData;
    tags = Set(union(tags, [:FACTOR])),
    solvable,
    timestamp = _zonedtime(timestamp),
  )
  #

  factor = addFactor!(dfg, newFactor)

  # TODO: change this operation to update a conditioning variable
  graphinit && doautoinit!(dfg, Xi; singles = false)

  return factor
end

function _checkFactorAdd(usrfnc, xisyms)
  if length(xisyms) == 1 && !(usrfnc isa AbstractPrior) && !(usrfnc isa Mixture)
    @warn("Listing only one variable $xisyms for non-unary factor type $(typeof(usrfnc))")
  end
  return nothing
end

function DFG.addFactor!(
  dfg::AbstractDFG,
  vlbs::AbstractVector{Symbol},
  usrfnc::AbstractFactor;
  suppressChecks::Bool = false,
  kw...,
)
  #

  # basic sanity check for unary vs n-ary
  if !suppressChecks
    _checkFactorAdd(usrfnc, vlbs)
    @assert length(vlbs) == length(unique(vlbs)) "List of variables should be unique and ordered."
  end

  # variables = getVariable.(dfg, vlbs)
  variables = map(vid -> getVariable(dfg, vid), vlbs)
  return addFactor!(dfg, variables, usrfnc; suppressChecks, kw...)
end

#
