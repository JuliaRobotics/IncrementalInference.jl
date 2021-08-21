


"""
$SIGNATURES

Initialize an empty in-memory DistributedFactorGraph `::DistributedFactorGraph` object.
"""
function initfg(dfg::T=InMemDFGType(solverParams=SolverParams());
                                sessionname="NA",
                                robotname="",
                                username="",
                                cloudgraph=nothing) where T <: AbstractDFG
#
  #
  return dfg
end


#init an empty fg with a provided type and SolverParams
function initfg(::Type{T};solverParams=SolverParams(),
                      sessionname="NA",
                      robotname="",
                      username="",
                      cloudgraph=nothing) where T <: AbstractDFG
  #
  return T(solverParams=solverParams)
end

function initfg(::Type{T},solverParams::S;
                      sessionname="NA",
                      robotname="",
                      username="",
                      cloudgraph=nothing) where {T <: AbstractDFG, S <: SolverParams}
  #
  return T{S}(solverParams=solverParams)
end



# Should deprecate in favor of TensorCast.jl
reshapeVec2Mat(vec::Vector, rows::Int) = reshape(vec, rows, round(Int,length(vec)/rows))



## ==============================================================================================
## MOVE TO / CONSOLIDATE WITH DFG
## ==============================================================================================


"""
    $(SIGNATURES)

Fetch the variable marginal joint sampled points.  Use [`getBelief`](@ref) to retrieve the full Belief object.
"""
getVal(v::DFGVariable; solveKey::Symbol=:default) = v.solverDataDict[solveKey].val
getVal(v::DFGVariable, idx::Int; solveKey::Symbol=:default) = v.solverDataDict[solveKey].val[:,idx]
getVal(vnd::VariableNodeData) = vnd.val
getVal(vnd::VariableNodeData, idx::Int) = vnd.val[:, idx]
getVal(dfg::AbstractDFG, lbl::Symbol; solveKey::Symbol=:default) = getVariable(dfg, lbl).solverDataDict[solveKey].val

"""
    $(SIGNATURES)

Get the number of points used for the current marginal belief estimate represtation for a particular variable in the factor graph.
"""
function getNumPts(v::DFGVariable; solveKey::Symbol=:default)::Int
  return length(getVal(getSolverData(v, solveKey)))
end

function AMP.getBW(vnd::VariableNodeData)
  return vnd.bw
end

# setVal! assumes you will update values to database separate, this used for local graph mods only
function getBWVal(v::DFGVariable; solveKey::Symbol=:default)
  return getSolverData(v, solveKey).bw
end
function setBW!(vd::VariableNodeData, bw::Array{Float64,2}; solveKey::Symbol=:default)
  vd.bw = bw
  nothing
end
function setBW!(v::DFGVariable, bw::Array{Float64,2}; solveKey::Symbol=:default)
  setBW!(getSolverData(v, solveKey), bw)
  nothing
end

function setVal!(vd::VariableNodeData, val::AbstractVector{P}) where P
  vd.val = val
  nothing
end
function setVal!(v::DFGVariable, val::AbstractVector{P}; solveKey::Symbol=:default) where P
    setVal!(getSolverData(v, solveKey), val)
    nothing
end
function setVal!(vd::VariableNodeData, val::AbstractVector{P}, bw::Array{Float64,2}) where P
    setVal!(vd, val)
    setBW!(vd, bw)
    nothing
end
function setVal!(v::DFGVariable, val::AbstractVector{P}, bw::Array{Float64,2}; solveKey::Symbol=:default) where P
  setVal!(v, val, solveKey=solveKey)
  setBW!(v, bw, solveKey=solveKey)
  nothing
end
function setVal!(vd::VariableNodeData, val::AbstractVector{P}, bw::Vector{Float64}) where P
  setVal!(vd, val, reshape(bw,length(bw),1))
  nothing
end
function setVal!(v::DFGVariable, val::AbstractVector{P}, bw::Vector{Float64}; solveKey::Symbol=:default) where P
  setVal!(getSolverData(v, solveKey), val, bw)
  nothing
end
function setVal!(dfg::AbstractDFG, sym::Symbol, val::AbstractVector{P}; solveKey::Symbol=:default) where P
  setVal!(getVariable(dfg, sym), val, solveKey=solveKey)
end

"""
    $SIGNATURES

Set the point centers and bandwidth parameters of a variable node, also set `isInitialized=true` if `setinit::Bool=true` (as per default).

Notes
- `initialized` is used for initial solve of factor graph where variables are not yet initialized.
- `inferdim` is used to identify if the initialized was only partial.
"""
function setValKDE!(vd::VariableNodeData,
                    pts::AbstractVector{P},
                    bws::Vector{Float64},
                    setinit::Bool=true,
                    inferdim::Float64=0.0 ) where P
  #

  
  setVal!(vd, pts, bws) # BUG ...al!(., val, . ) ## TODO -- this can be a little faster
  setinit ? (vd.initialized = true) : nothing
  vd.inferdim = inferdim
  nothing
end

function setValKDE!(vd::VariableNodeData,
                    val::AbstractVector{P},
                    setinit::Bool=true,
                    inferdim::Real=0.0  ) where P
  # recover variableType information
  varType = getVariableType(vd)
  p = AMP.manikde!(val, varType)
  setValKDE!(vd, p, setinit, inferdim)
  nothing
end

function setValKDE!(v::DFGVariable,
                    val::AbstractVector{P},
                    bws::Array{<:Real,2},
                    setinit::Bool=true,
                    inferdim::Float64=0;
                    solveKey::Symbol=:default) where P
  # recover variableType information
  setValKDE!(getSolverData(v, solveKey), val, bws[:,1], setinit, inferdim )

  nothing
end

function setValKDE!(v::DFGVariable,
                    val::AbstractVector{P},
                    setinit::Bool=true,
                    inferdim::Float64=0.0;
                    solveKey::Symbol=:default) where P
  # recover variableType information
  setValKDE!(getSolverData(v, solveKey),val, setinit, inferdim )
  nothing
end
function setValKDE!(v::DFGVariable,
                    em::TreeBelief,
                    setinit::Bool=true;
                    # inferdim::Union{Float32, Float64, Int32, Int64}=0;
                    solveKey::Symbol=:default  )
  #
  setValKDE!(v, em.val, em.bw, setinit, em.inferdim, solveKey=solveKey)
  nothing
end
function setValKDE!(v::DFGVariable,
                    mkd::ManifoldKernelDensity,
                    setinit::Bool=true,
                    inferdim::Real=0;
                    solveKey::Symbol=:default  )
  #
  # @error("TESTING setValKDE! ", solveKey, string(listSolveKeys(v)))
  setValKDE!(getSolverData(v,solveKey),mkd,setinit, Float64(inferdim))
  nothing
end
function setValKDE!(dfg::AbstractDFG,
                    sym::Symbol,
                    mkd::ManifoldKernelDensity,
                    setinit::Bool=true,
                    inferdim::Real=0;
                    solveKey::Symbol=:default  )
  #
  setValKDE!(getVariable(dfg, sym), mkd, setinit, inferdim, solveKey=solveKey)
  nothing
end



function setValKDE!(vnd::VariableNodeData,
                    mkd::ManifoldKernelDensity{M,B,Nothing},
                    setinit::Bool=true,
                    inferdim::Union{Float32, Float64, Int32, Int64}=0 ) where {M,B}
  #
  # L==Nothing means no partials
  ptsArr = AMP.getPoints(mkd) # , false) # for not partial
  # also set the bandwidth
  bws = getBW(mkd)[:,1]
  setValKDE!(vnd,ptsArr,bws,setinit,inferdim )
  nothing
end


function setValKDE!(vnd::VariableNodeData,
                    mkd::ManifoldKernelDensity{M,B,L},
                    setinit::Bool=true,
                    inferdim::Union{Float32, Float64, Int32, Int64}=0 ) where {M,B,L<:AbstractVector}
  #
  oldBel = getBelief(vnd)

  # New infomation might be partial
  newBel = replace(oldBel, mkd)

  # Set partial dims as Manifold points
  ptsArr = AMP.getPoints(newBel, false)

  # also get the bandwidth
  bws = getBandwidth(newBel, false)

  # update values in graph
  setValKDE!(vnd,ptsArr,bws,setinit,inferdim )
  nothing
end


"""
    $SIGNATURES

Set variable initialized status.
"""
function setVariableInitialized!( varid::VariableNodeData,
                                  status::Bool)
  #
  varid.initialized = status
end

setVariableInitialized!(vari::DFGVariable, status::Bool) = setVariableInitialized!(getSolverData(vari), status)


"""
    $SIGNATURES

Set method for the inferred dimension value in a variable.
"""
setVariableInferDim!(varid::VariableNodeData, val::Real) = varid.inferdim = convert(Float64,val)
setVariableInferDim!(vari::DFGVariable, val::Real) = setVariableInferDim!(getSolverData(vari), val)

## ==============================================================================================
## ==============================================================================================

"""
    $(SIGNATURES)

Get a ManifoldKernelDensity estimate from variable node data.
"""
getBelief(vnd::VariableNodeData) = manikde!(getManifold(getVariableType(vnd)), getVal(vnd), bw=getBW(vnd)[:,1] )

getBelief(v::DFGVariable, solvekey::Symbol=:default) = getBelief(getSolverData(v, solvekey))
getBelief(dfg::AbstractDFG, lbl::Symbol, solvekey::Symbol=:default) = getBelief(getVariable(dfg, lbl), solvekey)


"""
    $SIGNATURES

Reset the solve state of a variable to uninitialized/unsolved state.
"""
function resetVariable!(varid::VariableNodeData;
                        solveKey::Symbol=:default  )::Nothing
  #
  val = getBelief(varid)
  pts = AMP.getPoints(val)
  # TODO not all manifolds will initialize to zero
  for pt in pts
    fill!(pt, 0.0)
  end
  pn = manikde!(getManifold(varid), pts, bw=zeros(Ndim(val)))
  setValKDE!(varid, pn, false, 0.0)
  # setVariableInferDim!(varid, 0)
  # setVariableInitialized!(vari, false)
  nothing
end

resetVariable!(vari::DFGVariable; solveKey::Symbol=:default  ) = resetVariable!(getSolverData(vari), solveKey=solveKey)

function resetVariable!(dfg::G,
                        sym::Symbol;
                        solveKey::Symbol=:default  )::Nothing where G <: AbstractDFG
  #
  resetVariable!(getVariable(dfg, sym), solveKey=solveKey)
end




# TODO -- there should be a better way, without retrieving full vertex
# TODO -- Deprecated for DFG -- must update
function getOutNeighbors(dfg::T, v::V; needdata::Bool=false, ready::Union{Nothing, Int}=nothing, backendset::Union{Nothing, Int}=nothing)::Vector{Symbol} where {T <: AbstractDFG, V <: DFGNode}
  @warn "TODO: needdata is currently ignored. Symbols are returned."
  nodes = getNeighbors(dfg, v, ready=ready, backendset=backendset)
  return nodes
end
function getOutNeighbors(dfg::T, vertSym::Symbol; needdata::Bool=false, ready::Int=1, backendset::Int=1 )::Vector{Symbol} where {T <: AbstractDFG, V <: DFGNode}
  @warn "TODO: needdata is currently ignored. Symbols are returned."
  nodes = getNeighbors(dfg, vertSym, ready=ready, backendset=backendset)
  return nodes
end



function DefaultNodeDataParametric( dodims::Int,
                                    dims::Int,
                                    variableType::InferenceVariable;
                                    initialized::Bool=true,
                                    dontmargin::Bool=false)::VariableNodeData

  # this should be the only function allocating memory for the node points
  if false && initialized
    error("not implemented yet")
    # pN = AMP.manikde!(randn(dims, N), variableType.manifolds);
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
    sp = round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
    ϵ = getPointIdentity(variableType)
    return VariableNodeData([ϵ],
                            zeros(dims,dims), Symbol[], sp,
                            dims, false, :_null, Symbol[], variableType, false, 0.0, false, dontmargin, 0, 0, :parametric)
  end

end

function setDefaultNodeDataParametric!(v::DFGVariable, variableType::InferenceVariable; kwargs...)
  vnd = DefaultNodeDataParametric(0, variableType |> getDimension, variableType; kwargs...)
  setSolverData!(v, vnd, :parametric)
  return nothing
end

"""
    $SIGNATURES

Create new solverData.

Notes
- Used during creation of new variable, as well as in CSM unique `solveKey`.
"""
function setDefaultNodeData!( v::DFGVariable,
                              dodims::Int,
                              N::Int,
                              dims::Int;
                              solveKey::Symbol=:default,
                              gt=Dict(),
                              initialized::Bool=true,
                              dontmargin::Bool=false,
                              varType=nothing )
  #
  # TODO review and refactor this function, exists as legacy from pre-v0.3.0
  # this should be the only function allocating memory for the node points (unless number of points are changed)
  data = nothing
  isinit = false
  sp = Int[0;]
  (valpts, bws) = if initialized
    pN = resample(getBelief(v))
    bws = getBW(pN)[:,1:1]
    pNpts = getPoints(pN)
    isinit = true
    (pNpts, bws)
  else
    sp = round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
    @assert getPointType(varType) != DataType "cannot add manifold point type $(getPointType(varType)), make sure the identity element argument in @defVariable $varType arguments is correct"
    valpts = Vector{getPointType(varType)}(undef,N)
    for i in 1:length(valpts)
      valpts[i] = getPointIdentity(varType)
    end
    bws = zeros(dims,1)
    #
    (valpts, bws)
  end
  # make and set the new solverData
  setSolverData!(v, VariableNodeData( valpts, bws,
                                      Symbol[], sp,
                                      dims, false, :_null, Symbol[], 
                                      varType, isinit, 0.0, false, dontmargin,0,0,solveKey), solveKey)
  return nothing
end
# if size(initval,2) < N && size(initval, 1) == dims
#   @warn "setDefaultNodeData! -- deprecated use of stdev."
#   p = AMP.manikde!(initval,diag(stdev), varType.manifolds);
#   pN = resample(p,N)
# if size(initval,2) < N && size(initval, 1) != dims
  # @info "Node value memory allocated but not initialized"
# else
#   pN = AMP.manikde!(initval, varType.manifolds)
# end
# dims = size(initval,1) # rows indicate dimensions



"""
    $SIGNATURES

Reference data can be stored in the factor graph as a super-solve.

Notes
- Intended as a mechanism to store reference data alongside the numerical computations.
"""
function setVariableRefence!( dfg::AbstractDFG,
                              sym::Symbol,
                              val::AbstractVector;
                              refKey::Symbol=:reference)
  #
  # which variable to update
  var = getVariable(dfg, sym)

  # Construct an empty VND object
  vnd = VariableNodeData( val,
                          zeros(getDimension(var),1),
                          Symbol[],
                          Int[0;],
                          getDimension(var),
                          false,
                          :_null,
                          Symbol[],
                          getVariableType(var),
                          true,
                          0.0,
                          false,
                          true  )
  #
  # set the value in the DFGVariable
  setSolverData!(var, vnd, refKey)
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
function addVariable!(dfg::AbstractDFG,
                      label::Symbol,
                      varTypeU::Union{T, Type{T}}; 
                      N::Int=getSolverParams(dfg).N,
                      solvable::Int=1,
                      timestamp::Union{DateTime,ZonedDateTime}=now(localzone()),
                      nanosecondtime::Union{Nanosecond,Int64,Nothing}=nothing,
                      dontmargin::Bool=false,
                      labels::Union{Vector{Symbol},Nothing}=nothing,
                      tags::Vector{Symbol}=Symbol[],
                      smalldata=Dict{Symbol, DFG.SmallDataTypes}(),
                      checkduplicates::Bool=true,
                      initsolvekeys::Vector{Symbol}=getSolverParams(dfg).algorithms ) where T<:InferenceVariable
  #
  varType = _variableType(varTypeU)
  # TODO Remove deprecation in v0.16 
  if :ut in fieldnames(T)	
    Base.depwarn("Field `ut` (microseconds) for variable type ($T) has been deprecated please use DFGVariable.nstime, kwarg: nanosecondtime", :addVariable!)	
    if isnothing(nanosecondtime)
      varType.ut == -9999999999 && error("please define a time for type $(T), use FGVariable.nstime, kwarg: nanosecondtime")
      nanosecondtime = Nanosecond(varType.ut*1000)	
    else 	
      @warn "Nanosecond time has been specified as $nanosecondtime, ignoring `ut` field value: $(varType.ut)."	
    end	
  elseif isnothing(nanosecondtime)	
    nanosecondtime = Nanosecond(0)	
  end

  # deprecate in v0.16
  labels isa Vector ? (union!(tags, labels); @warn("labels is deprecated, use tags instead")) : nothing
  
  union!(tags, [:VARIABLE])
  v = DFGVariable(label, varType; tags=Set(tags), smallData=smalldata, solvable=solvable, timestamp=timestamp, nstime=Nanosecond(nanosecondtime))

  (:default in initsolvekeys) &&
    setDefaultNodeData!(v, 0, N, getDimension(varType), initialized=false, varType=varType, dontmargin=dontmargin) # dodims

  (:parametric in initsolvekeys) &&
    setDefaultNodeDataParametric!(v, varType, initialized=false, dontmargin=dontmargin)

  DFG.addVariable!(dfg, v)

  return v
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
    multihypo2[1-1e-10 .< multihypo] .= 0.0
    # check that terms sum to full probability
    @assert abs(sum(multihypo2) % 1) < 1e-10  || 1-1e-10 < sum(multihypo2) % 1 "ensure multihypo sums to a (or nearly, 1e-10) interger, see #1086"
    # check that only one variable broken into fractions
    @assert sum(multihypo2[1e-10 .< multihypo2]) ≈ 1
    # force normalize something that is now known to be close
    multihypo2 ./= sum(multihypo2)

    mh = Categorical(Float64[multihypo2...] )
  end
  return mh, nullhypo
end



# return a BitVector masking the fractional portion, assuming converted 0's on 100% confident variables 
_getFractionalVars(varList::Union{<:Tuple, <:AbstractVector}, mh::Nothing) = zeros(length(varList)) .== 1
_getFractionalVars(varList::Union{<:Tuple, <:AbstractVector}, mh::Categorical) = 0 .< mh.p

function _selectHypoVariables(allVars::Union{<:Tuple, <:AbstractVector}, 
                              mh::Categorical,
                              sel::Integer = rand(mh) )
  #
  mask = mh.p .≈ 0.0
  mask[sel] = true
  (1:length(allVars))[mask]
end

_selectHypoVariables(allVars::Union{<:Tuple, <:AbstractVector},mh::Nothing,sel::Integer=0 ) = collect(1:length(allVars))


function prepgenericconvolution(Xi::Vector{<:DFGVariable},
                                usrfnc::T;
                                multihypo::Union{Nothing, Distributions.Categorical}=nothing,
                                nullhypo::Real=0.0,
                                threadmodel=MultiThreaded,
                                inflation::Real=0.0,
                                _blockRecursion::Bool=false  ) where {T <: AbstractFactor}
  #
  pttypes = getVariableType.(Xi) .|> getPointType
  PointType = 0 < length(pttypes) ? pttypes[1] : Vector{Float64}
  # FIXME stop using Any, see #1321
  varParamsAll = Vector{Vector{Any}}()
  maxlen, sfidx, mani = prepareparamsarray!(varParamsAll, Xi, nothing, 0) # Nothing for init.

  # standard factor metadata
  sflbl = 0==length(Xi) ? :null : getLabel(Xi[end])
  fmd = FactorMetadata(Xi, getLabel.(Xi), varParamsAll, sflbl, nothing)
  
  # create a temporary CalcFactor object for extracting the first sample
  # TODO, deprecate this:  guess measurement points type
  # MeasType = Vector{Float64} # FIXME use `usrfnc` to get this information instead
  _cf = CalcFactor( usrfnc, fmd, 0, 1, nothing, varParamsAll) # (Vector{MeasType}(),)
  
  # get a measurement sample
  meas_single = sampleFactor(_cf, 1)

  #TODO preallocate measuerement?
  measurement = Vector{eltype(meas_single)}()
  
  # get the measurement dimension
  zdim = calcZDim(_cf)
  # some hypo resolution
  certainhypo = multihypo !== nothing ? collect(1:length(multihypo.p))[multihypo.p .== 0.0] : collect(1:length(Xi))
  
  # sort out partialDims here
  ispartl = hasfield(T, :partial)
  partialDims = if ispartl
    Int[usrfnc.partial...]
  else
    Int[]
  end

  # as per struct CommonConvWrapper
  varTypes::Vector{DataType} = typeof.(getVariableType.(Xi))
  gradients = nothing
  # prepare new cached gradient lambdas (attempt)
  try
    # https://github.com/JuliaRobotics/IncrementalInference.jl/blob/db7ff84225cc848c325e57b5fb9d0d85cb6c79b8/src/DispatchPackedConversions.jl#L46
    # also https://github.com/JuliaRobotics/DistributedFactorGraphs.jl/issues/590#issuecomment-891450762
    # FIXME, suppressing nested gradient propagation on GenericMarginals for the time being, see #1010
    if (!_blockRecursion) && usrfnc isa AbstractRelative && !(usrfnc isa GenericMarginal)
      # take first value from each measurement-tuple-element
      measurement_ = map(x->x[1], meas_single)
      # compensate if no info available during deserialization
      # take the first value from each variable param
      pts_ = map(x->x[1], varParamsAll)
      # FIXME, only using first meas and params values at this time...
      # NOTE, must block recurions here, since FGC uses this function to calculate numerical gradients on a temp fg.
      # assume for now fractional-var in multihypo have same varType
      hypoidxs = _selectHypoVariables(pts_, multihypo)
      gradients = FactorGradientsCached!(usrfnc, tuple(varTypes[hypoidxs]...), measurement_, tuple(pts_[hypoidxs]...), _blockRecursion=true);
    end
  catch e
    @warn "Unable to create measurements and gradients for $usrfnc during prep of CCW, falling back on no-partial information assumption.  Enable ENV[\"JULIA_DEBUG\"] = \"IncrementalInference\" for @debug printing to see the error."
    @debug(e)
  end

  ccw = CommonConvWrapper(
          usrfnc,
          PointType[],
          zdim,
          varParamsAll,
          fmd;
          specialzDim = hasfield(T, :zDim),
          partial = ispartl,
          measurement,
          hypotheses=multihypo,
          certainhypo=certainhypo,
          nullhypo=nullhypo,
          threadmodel=threadmodel,
          inflation=inflation,
          partialDims=partialDims,
          vartypes = varTypes,
          gradients=gradients
        )
  #
  return ccw
end

# TODO perhaps consolidate with constructor?
"""
$SIGNATURES

Generate the default factor data for a new DFGFactor.
"""
function getDefaultFactorData(dfg::AbstractDFG,
                              Xi::Vector{<:DFGVariable},
                              usrfnc::T;
                              multihypo::Vector{<:Real}=Float64[],
                              nullhypo::Float64=0.0,
                              threadmodel=SingleThreaded,
                              eliminated::Bool = false,
                              potentialused::Bool = false,
                              edgeIDs = Int[],
                              solveInProgress = 0,
                              inflation::Real=getSolverParams(dfg).inflation,
                              _blockRecursion::Bool=false ) where T <: FunctorInferenceType
  #

  # prepare multihypo particulars
  # storeMH::Vector{Float64} = multihypo == nothing ? Float64[] : [multihypo...]
  mhcat, nh = parseusermultihypo(multihypo, nullhypo)

  # allocate temporary state for convolutional operations (not stored)
  ccw = prepgenericconvolution(Xi, usrfnc, multihypo=mhcat, nullhypo=nh, threadmodel=threadmodel, inflation=inflation, _blockRecursion=_blockRecursion)

  # and the factor data itself
  return FunctionNodeData{typeof(ccw)}(eliminated, potentialused, edgeIDs, ccw, multihypo, ccw.certainhypo, nullhypo, solveInProgress, inflation)
end


"""
    $SIGNATURES

Return `::Bool` on whether at least one hypothesis is available for intended computations (assuming direction `sfidx`).
"""
function isLeastOneHypoAvailable( sfidx::Int,
                                  certainidx::Vector{Int},
                                  uncertnidx::Vector{Int},
                                  isinit::Vector{Bool})
  #
  # @show isinit
  # @show sfidx in certainidx, sum(isinit[uncertnidx])
  # @show sfidx in uncertnidx, sum(isinit[certainidx])
  return  sfidx in certainidx && 0 < sum(isinit[uncertnidx]) ||
          sfidx in uncertnidx && sum(isinit[certainidx]) == length(certainidx)
end


function assembleFactorName(dfg::AbstractDFG,
                            Xi::Vector{<:DFGVariable} )
  #

  existingFactorLabels = listFactors(dfg)
  existingFactorLabelDict = Dict(existingFactorLabels .=> existingFactorLabels)
  namestring = ""
  for vert in Xi #f.Xi
    namestring = string(namestring,vert.label)
  end
  opt = getSolverParams(dfg)
  for i in 1:opt.maxincidence
    tempnm = string(namestring, "f$i")
    if !haskey(existingFactorLabelDict, Symbol(tempnm))
      namestring = tempnm
      break
    end
    i != opt.maxincidence ? nothing : error("Artificial restriction to not connect more than $(opt.maxincidence) factors to a variable (bad for sparsity), try setting getSolverParams(fg).maxincidence=1000 to adjust this restriction.")
  end
  return Symbol(namestring)
end



"""
    $(SIGNATURES)

Add factor with user defined type <: FunctorInferenceType to the factor graph
object. Define whether the automatic initialization of variables should be
performed.  Use order sensitive `multihypo` keyword argument to define if any
variables are related to data association uncertainty.

Experimental
- `inflation`, to better disperse kernels before convolution solve, see IIF #1051.
"""
function DFG.addFactor!(dfg::AbstractDFG,
                        Xi::AbstractVector{<:DFGVariable},
                        usrfnc::AbstractFactor;
                        multihypo::Vector{Float64}=Float64[],
                        nullhypo::Float64=0.0,
                        solvable::Int=1,
                        tags::Vector{Symbol}=Symbol[],
                        timestamp::Union{DateTime,ZonedDateTime}=now(localzone()),
                        graphinit::Bool=getSolverParams(dfg).graphinit,
                        threadmodel=SingleThreaded,
                        suppressChecks::Bool=false,
                        inflation::Real=getSolverParams(dfg).inflation,
                        namestring::Symbol = assembleFactorName(dfg, Xi),
                        _blockRecursion::Bool=!getSolverParams(dfg).attemptGradients  )
  #

  varOrderLabels = Symbol[v.label for v=Xi]
  solverData = getDefaultFactorData(dfg, 
                                    Xi, 
                                    deepcopy(usrfnc), 
                                    multihypo=multihypo, 
                                    nullhypo=nullhypo, 
                                    threadmodel=threadmodel,
                                    inflation=inflation,
                                    _blockRecursion=_blockRecursion)
  newFactor = DFGFactor(Symbol(namestring),
                        varOrderLabels,
                        solverData;
                        tags=Set(union(tags, [:FACTOR])),
                        solvable=solvable,
                        timestamp=timestamp)
  #

  success = addFactor!(dfg, newFactor)

  # TODO: change this operation to update a conditioning variable
  graphinit && doautoinit!(dfg, Xi, singles=false)

  return newFactor
end

function _checkFactorAdd(usrfnc, xisyms)
  if length(xisyms) == 1 && !(usrfnc isa AbstractPrior) && !(usrfnc isa Mixture)
    @warn("Listing only one variable $xisyms for non-unary factor type $(typeof(usrfnc))")
  end
  nothing
end

function DFG.addFactor!(dfg::AbstractDFG,
                        xisyms::AbstractVector{Symbol},
                        usrfnc::AbstractFactor;
                        suppressChecks::Bool=false,
                        kw...  )
  #

  # basic sanity check for unary vs n-ary
  if !suppressChecks
    _checkFactorAdd(usrfnc, xisyms)
  end

  # variables = getVariable.(dfg, xisyms)
  variables = map(vid -> getVariable(dfg, vid), xisyms)
  addFactor!(dfg, variables, usrfnc; suppressChecks=suppressChecks, kw... ) 
end





#
