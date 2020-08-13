
# Should deprecate in favor of TensorCast.jl
reshapeVec2Mat(vec::Vector, rows::Int) = reshape(vec, rows, round(Int,length(vec)/rows))

"""
    $SIGNATURES
Return the manifolds on which variable `sym::Symbol` is defined.
"""
getManifolds(vd::VariableNodeData) = getSofttype(vd) |> getManifolds
getManifolds(v::DFGVariable; solveKey::Symbol=:default) = getManifolds(getSolverData(v, solveKey))
function getManifolds(dfg::G, sym::Symbol; solveKey::Symbol=:default) where {G <: AbstractDFG}
  return getManifolds(getVariable(dfg, sym), solveKey=solveKey)
end
getManifolds(vartype::InferenceVariable) = vartype.manifolds
getManifolds(vartype::Type{<: InferenceVariable}) = getManifolds(vartype())

"""
    $(SIGNATURES)

Convenience function to get point values sampled i.i.d from marginal of `lbl` variable in the current factor graph.
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
  return size(getSolverData(v, solveKey).val,2)
end

# import DistributedFactorGraphs: getfnctype
# # TODO: Refactor - was is das?
# function getfnctype(data::GenericFunctionNodeData)
#   if typeof(data).name.name == :VariableNodeData
#     return VariableNodeData
#   end
#   return data.fnc.usrfnc!
# end
#
# function getfnctype(fact::DFGFactor; solveKey::Symbol=:default)
#   data = getData(fact) # TODO , solveKey=solveKey)
#   return getfnctype(data)
# end
#
# function getfnctype(dfg::T, lbl::Symbol; solveKey::Symbol=:default) where T <: AbstractDFG
#   getfnctype(getFactor(dfg, exvertid))
# end

function getBW(vnd::VariableNodeData)
  return vnd.bw
end

# setVal! assumes you will update values to database separate, this used for local graph mods only
function getBWVal(v::DFGVariable; solveKey::Symbol=:default)
  return getSolverData(v, solveKey).bw
end
function setBW!(vd::VariableNodeData, bw::Array{Float64,2}; solveKey::Symbol=:default)::Nothing
  vd.bw = bw
  nothing
end
function setBW!(v::DFGVariable, bw::Array{Float64,2}; solveKey::Symbol=:default)::Nothing
  setBW!(getSolverData(v, solveKey), bw)
  nothing
end

function setVal!(vd::VariableNodeData, val::Array{Float64,2})::Nothing
    vd.val = val
    nothing
end
function setVal!(v::DFGVariable, val::Array{Float64,2}; solveKey::Symbol=:default)::Nothing
    setVal!(getSolverData(v, solveKey), val)
    nothing
end
function setVal!(vd::VariableNodeData, val::Array{Float64,2}, bw::Array{Float64,2})::Nothing
    setVal!(vd, val)
    setBW!(vd, bw)
    nothing
end
function setVal!(v::DFGVariable, val::Array{Float64,2}, bw::Array{Float64,2}; solveKey::Symbol=:default)::Nothing
  setVal!(v, val, solveKey=solveKey)
  setBW!(v, bw, solveKey=solveKey)
  nothing
end
function setVal!(vd::VariableNodeData, val::Array{Float64,2}, bw::Vector{Float64}; solveKey::Symbol=:default)
  setVal!(vd, val, reshape(bw,length(bw),1))
  nothing
end
function setVal!(v::DFGVariable, val::Array{Float64,2}, bw::Vector{Float64}; solveKey::Symbol=:default)
  setVal!(getSolverData(v, solveKey=solveKey), val, bw)
  nothing
end
function setVal!(dfg::AbstractDFG, sym::Symbol, val::Array{Float64,2}; solveKey::Symbol=:default)
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
                    pts::Array{Float64,2},
                    bws::Vector{Float64},
                    setinit::Bool=true,
                    inferdim::Float64=0.0)::Nothing
  #
  setVal!(vd, pts, bws) # BUG ...al!(., val, . ) ## TODO -- this can be a little faster
  setinit ? (vd.initialized = true) : nothing
  vd.inferdim = inferdim
  nothing
end

function setValKDE!(vd::VariableNodeData,
                    p::BallTreeDensity,
                    setinit::Bool=true,
                    inferdim::Union{Float32, Float64, Int32, Int64}=0 )
  #
  pts = getPoints(p)
  bws = getBW(p)[:,1]
  setValKDE!(vd,pts,bws,setinit,inferdim )
  nothing
end

function setValKDE!(vd::VariableNodeData,
                    val::Array{Float64,2},
                    setinit::Bool=true,
                    inferdim::Float64=0.0)::Nothing
  # recover softtype information
  sty = getSofttype(vd)
  p = AMP.manikde!(val, getManifolds(sty))
  setValKDE!(vd, p, setinit, inferdim)
  nothing
end

function setValKDE!(v::DFGVariable,
                    val::Array{Float64,2},
                    bws::Array{Float64,2},
                    setinit::Bool=true,
                    inferdim::Float64=0;
                    solveKey::Symbol=:default)::Nothing
  # recover softtype information
  setValKDE!(getSolverData(v, solveKey), val, bws[:,1], setinit, inferdim )

  nothing
end

function setValKDE!(v::DFGVariable,
                    val::Array{Float64,2},
                    setinit::Bool=true,
                    inferdim::Float64=0.0;
                    solveKey::Symbol=:default)::Nothing
  # recover softtype information
  setValKDE!(getSolverData(v, solveKey),val, setinit, inferdim )
  nothing
end
function setValKDE!(v::DFGVariable,
                    em::TreeBelief,
                    setinit::Bool=true;
                    # inferdim::Union{Float32, Float64, Int32, Int64}=0;
                    solveKey::Symbol=:default  )::Nothing
  #
  setValKDE!(v, em.val, em.bw, setinit, em.inferdim, solveKey=solveKey)
  nothing
end
function setValKDE!(v::DFGVariable,
                    p::BallTreeDensity,
                    setinit::Bool=true,
                    inferdim::Union{Float32, Float64, Int32, Int64}=0;
                    solveKey::Symbol=:default  )
  #
  setValKDE!(getSolverData(v,solveKey),p,setinit,Float64(inferdim))
  nothing
end
function setValKDE!(dfg::G,
                    sym::Symbol,
                    p::BallTreeDensity,
                    setinit::Bool=true,
                    inferdim::Union{Float32, Float64, Int32, Int64}=0;
                    solveKey::Symbol=:default  ) where G <: AbstractDFG
    #
    setValKDE!(getVariable(dfg, sym), p, setinit, inferdim, solveKey=solveKey)
    nothing
end


"""
    $SIGNATURES

Set variable initialized status.
"""
function setVariableInitialized!(varid::VariableNodeData,
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


"""
    $SIGNATURES

Reset the solve state of a variable to uninitialized/unsolved state.
"""
function resetVariable!(varid::VariableNodeData;
                        solveKey::Symbol=:default  )::Nothing
  #
  val = getKDE(varid)
  pts = getPoints(val)
  # TODO not all manifolds will initialize to zero
  fill!(pts, 0.0)
  pn = manikde!(pts, zeros(KDE.Ndim(val)), getManifolds(varid))
  setValKDE!(varid, pn, false, 0.0)
  # setVariableInferDim!(varid, 0)
  # setVariableInitialized!(vari, false)
  nothing
end

resetVariable!(vari::DFGVariable; solveKey::Symbol=:default  )::Nothing = resetVariable!(getSolverData(vari), solveKey=solveKey)

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



function DefaultNodeDataParametric(dodims::Int,
                                   dims::Int,
                                   softtype::InferenceVariable;
                                   initialized::Bool=true,
                                   dontmargin::Bool=false)::VariableNodeData

  # this should be the only function allocating memory for the node points
  if false && initialized
    error("not implemented yet")
    # pN = AMP.manikde!(randn(dims, N), softtype.manifolds);
    #
    # sp = Int[0;] #round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
    # gbw = getBW(pN)[:,1]
    # gbw2 = Array{Float64}(undef, length(gbw),1)
    # gbw2[:,1] = gbw[:]
    # pNpts = getPoints(pN)
    # #initval, stdev
    # return VariableNodeData(pNpts,
    #                         gbw2, Symbol[], sp,
    #                         dims, false, :_null, Symbol[], softtype, true, 0.0, false, dontmargin)
  else
    sp = round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
    return VariableNodeData(zeros(dims, 1),
                            zeros(dims,1), Symbol[], sp,
                            dims, false, :_null, Symbol[], softtype, false, 0.0, false, dontmargin)
  end

end

function setDefaultNodeDataParametric!(v::DFGVariable, softtype::InferenceVariable; kwargs...)
  vnd = DefaultNodeDataParametric(0, softtype |> getDimension, softtype; kwargs...)
  setSolverData!(v, vnd, :parametric)
  return nothing
end

function setDefaultNodeData!(v::DFGVariable,
                             dodims::Int,
                             N::Int,
                             dims::Int;
                             gt=Dict(),
                             initialized::Bool=true,
                             dontmargin::Bool=false,
                             varType=nothing)::Nothing
  # TODO review and refactor this function, exists as legacy from pre-v0.3.0
  # this should be the only function allocating memory for the node points (unless number of points are changed)
  data = nothing
  if initialized

      pN = AMP.manikde!(randn(dims, N), getManifolds(varType));

    sp = Int[0;] #round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
    gbw = getBW(pN)[:,1]
    gbw2 = Array{Float64}(undef, length(gbw),1)
    gbw2[:,1] = gbw[:]
    pNpts = getPoints(pN)
    #initval, stdev
    setSolverData!(v, VariableNodeData(pNpts,
                            gbw2, Symbol[], sp,
                            dims, false, :_null, Symbol[], varType, true, 0.0, false, dontmargin))
  else
    sp = round.(Int,range(dodims,stop=dodims+dims-1,length=dims))
    setSolverData!(v, VariableNodeData(zeros(dims, N),
                            zeros(dims,1), Symbol[], sp,
                            dims, false, :_null, Symbol[], varType, false, 0.0, false, dontmargin))
  end
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
function setVariableRefence!(dfg::AbstractDFG,
                             sym::Symbol,
                             val::Array{Float64,2};
                             refKey::Symbol=:reference)
  #
  # which variable to update
  var = getVariable(dfg, sym)

  # Construct an empty VND object
  vnd = VariableNodeData(val,
                         zeros(getDimension(var),1),
                         Symbol[],
                         Int[0;],
                         getDimension(var),
                         false,
                         :_null,
                         Symbol[],
                         getSofttype(var),
                         true,
                         0.0,
                         false,
                         true  )
  #
  # set the value in the DFGVariable
  setSolverData!(var, vnd, refKey)
end


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
                      varType::InferenceVariable;
                      N::Int=100,
                      solvable::Int=1,
                      timestamp::Union{DateTime,ZonedDateTime}=now(localzone()),
                      nanosecondtime::Union{Nanosecond,Int64}=Nanosecond(0),
                      dontmargin::Bool=false,
                      labels::Union{Vector{Symbol},Nothing}=nothing,
                      tags::Vector{Symbol}=Symbol[],
                      smalldata=Dict{Symbol, DFG.SmallDataTypes}(),
                      checkduplicates::Bool=true,
                      initsolvekeys::Vector{Symbol}=getSolverParams(dfg).algorithms )
  #
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


function addVariable!(dfg::AbstractDFG,
                      label::Symbol,
                      stT::Type{<:InferenceVariable};
                      N::Int=100,
                      solvable::Int=1,
                      timestamp::Union{DateTime,ZonedDateTime}=now(localzone()),
                      nanosecondtime::Union{Nanosecond,Int64}=Nanosecond(0),
                      dontmargin::Bool=false,
                      labels::Union{Vector{Symbol},Nothing}=nothing,
                      tags::Vector{Symbol}=Symbol[],
                      smalldata=Dict{Symbol, DFG.SmallDataTypes}(),
                      checkduplicates::Bool=true,
                      initsolvekeys::Vector{Symbol}=getSolverParams(dfg).algorithms )
  #
  addVariable!(dfg,
               label,
               stT();
               N=N,
               solvable=solvable,
               timestamp=timestamp,
               nanosecondtime=nanosecondtime,
               dontmargin=dontmargin,
               labels=labels,
               tags=tags,
               smalldata=smalldata,
               checkduplicates=checkduplicates,
               initsolvekeys=initsolvekeys  )
end


"""
    $(SIGNATURES)

Fetch the variable marginal sample points without the KDE bandwidth parameter.  Use getVertKDE to retrieve the full KDE object.
"""
function getVal(vA::Vector{<:DFGVariable}, solveKey::Symbol=:default)::Array{Float64, 2}
  @warn "getVal(::Vector{DFGVariable}) is obsolete, use getVal.(DFGVariable) instead."
  len = length(vA)
  vals = Array{Array{Float64,2},1}()
  cols = Array{Int,1}()
  push!(cols,0)
  rows = Array{Int,1}()
  for v in vA
      push!(vals, getVal(v, solveKey=solveKey))
      c = size(vals[end],2)
      r = size(vals[end],1)
      push!(cols, floor(Int,c))
      push!(rows, floor(Int,r))
  end
  cols = cumsum(cols)
  sc = cols[end]
  rw = floor(Int,rows[1])
  val = Array{Float64,2}(undef,rw, sc)
  for i in 1:(len-1)
      val[:,(cols[i]+1):cols[i+1]] = vals[i]
  end
  val[:,(cols[len]+1):cols[len+1]] = vals[len] # and the last one
  return val
end


"""
    $(SIGNATURES)

Prepare the particle arrays `ARR` to be used for approximate convolution.
This function ensures that ARR has te same dimensions among all the parameters.
Function returns with ARR[sfidx] pointing at newly allocated deepcopy of the
existing values in getVal(Xi[.label==solvefor]).
Return values `sfidx` is the element in ARR where `Xi.label==solvefor` and
`maxlen` is length of all (possibly resampled) `ARR` contained particles.
Note `Xi` is order sensitive.
Note for initialization, solveFor = Nothing.
"""
function prepareparamsarray!(ARR::Array{Array{Float64,2},1},
                             Xi::Vector{<:DFGVariable},
                             solvefor::Union{Nothing, Symbol},
                             N::Int=0  )
  #
  LEN = Int[]
  maxlen = N # FIXME see #105
  count = 0
  sfidx = 0

  for xi in Xi
    push!(ARR, getVal(xi))
    len = size(ARR[end], 2)
    push!(LEN, len)
    if len > maxlen
      maxlen = len
    end
    count += 1
    if xi.label == solvefor
      sfidx = count #xi.index
    end
  end
  SAMP=LEN.<maxlen
  for i in 1:count
    if SAMP[i]
      ARR[i] = KDE.sample(getKDE(Xi[i]), maxlen)[1]
    end
  end

  # TODO --rather define reusable memory for the proposal
  # we are generating a proposal distribution, not direct replacement for existing memory and hence the deepcopy.
  if sfidx > 0 ARR[sfidx] = deepcopy(ARR[sfidx]) end

  # get solvefor manifolds
  manis = length(Xi)==0 || sfidx==0 ? (:null,) : getManifolds(Xi[sfidx])

  # FIXME, forcing maxlen to N results in errors (see test/testVariousNSolveSize.jl) see #105
  # maxlen = N == 0 ? maxlen : N
  return maxlen, sfidx, manis
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
    @assert sum(multihypo2) % 1 ≈ 0
    # check that only one variable broken into fractions
    @assert sum(multihypo2[1e-10 .< multihypo2]) ≈ 1
    # multihypo2 = Float64[multihypo...]
    # verts = Symbol.(multihypo[1,:])
    # for i in 1:length(multihypo)
    #   if multihypo[i] > 0.999999
    #     multihypo2[i] = 0.0
    #   end
    # end
    mh = Categorical(Float64[multihypo2...] )
  end
  return mh, nullhypo
end

# import IncrementalInference: prepgenericconvolution, convert

function calcZDim(usrfnc::T, Xi::Vector{<:DFGVariable})::Int where {T <: FunctorInferenceType}
  # zdim = T != GenericMarginal ? size(getSample(usrfnc, 2)[1],1) : 0
  zdim = if T != GenericMarginal
    vnds = Xi # (x->getSolverData(x)).(Xi)
    smpls = freshSamples(usrfnc, 2, FactorMetadata(), vnds)[1]
    size(smpls,1)
  else
    0
  end
  return zdim
end

function prepgenericconvolution(
            Xi::Vector{<:DFGVariable},
            usrfnc::T;
            multihypo::Union{Nothing, Distributions.Categorical}=nothing,
            nullhypo::Real=0.0,
            threadmodel=MultiThreaded  ) where {T <: FunctorInferenceType}
  #
  ARR = Array{Array{Float64,2},1}()
  maxlen, sfidx, manis = prepareparamsarray!(ARR, Xi, nothing, 0) # Nothing for init.
  fldnms = fieldnames(T) # typeof(usrfnc)
  zdim = calcZDim(usrfnc, Xi)
  # zdim = T != GenericMarginal ? size(getSample(usrfnc, 2)[1],1) : 0
  certainhypo = multihypo != nothing ? collect(1:length(multihypo.p))[multihypo.p .== 0.0] : collect(1:length(Xi))
  ccw = CommonConvWrapper(
          usrfnc,
          zeros(1,0),
          zdim,
          ARR,
          specialzDim = sum(fldnms .== :zDim) >= 1,
          partial = sum(fldnms .== :partial) >= 1,
          hypotheses=multihypo,
          certainhypo=certainhypo,
          nullhypo=nullhypo,
          threadmodel=threadmodel
        )
  #
  for i in 1:Threads.nthreads()
    # TODO JT - Confirm it should be updated here. Also testing in prepareCommonConvWrapper!
    ccw.cpt[i].factormetadata.fullvariables = copy(Xi)
    ccw.cpt[i].factormetadata.variableuserdata = []
    ccw.cpt[i].factormetadata.solvefor = :null
    for xi in Xi
      push!(ccw.cpt[i].factormetadata.variableuserdata, getSolverData(xi).softtype)
    end
  end
  return ccw
end

"""
$SIGNATURES

Generate the default factor data for a new DFGFactor.
"""
function getDefaultFactorData(
      dfg::AbstractDFG,
      Xi::Vector{<:DFGVariable},
      usrfnc::T;
      multihypo::Vector{<:Real}=Float64[],
      nullhypo::Float64=0.0,
      threadmodel=SingleThreaded ) where
        {T <: FunctorInferenceType}
  #

  # prepare multihypo particulars
  # storeMH::Vector{Float64} = multihypo == nothing ? Float64[] : [multihypo...]
  mhcat, nh = parseusermultihypo(multihypo, nullhypo)

  # allocate temporary state for convolutional operations (not stored)
  ccw = prepgenericconvolution(Xi, usrfnc, multihypo=mhcat, nullhypo=nh, threadmodel=threadmodel)

  # and the factor data itself
  data_ccw = FunctionNodeData{CommonConvWrapper{T}}(false, false, Int[], ccw, multihypo, ccw.certainhypo, nullhypo, 0)
  return data_ccw
end

"""
    $SIGNATURES

Returns state of vertex data `.initialized` flag.

Notes:
- used by Bayes tree clique logic.
- similar method in DFG
"""
function isInitialized(vert::TreeClique)::Bool
  return getSolverData(vert).initialized
end


"""
    $SIGNATURES

Return `::Bool` on whether at least one hypothesis is available for intended computations (assuming direction `sfidx`).
"""
function isLeastOneHypoAvailable(sfidx::Int,
                                 certainidx::Vector{Int},
                                 uncertnidx::Vector{Int},
                                 isinit::Vector{Bool})::Bool
  #
  # @show isinit
  # @show sfidx in certainidx, sum(isinit[uncertnidx])
  # @show sfidx in uncertnidx, sum(isinit[certainidx])
  return sfidx in certainidx && 0 < sum(isinit[uncertnidx]) ||
         sfidx in uncertnidx && sum(isinit[certainidx]) == length(certainidx)
end

"""
    $SIGNATURES

Return `(::Bool, ::OKVarlist, ::NotOkayVarList)` on whether all other variables (besides `loovar::Symbol`)
attached to factor `fct::Symbol` are all initialized -- i.e. `fct` is usable.

Notes:
- Special carve out for multihypo cases, see issue 427, where at least one hypothesis should be available, but not all required at first.

Development Notes
* TODO get faster version of isInitialized for database version

Related

doautoinit!, initManual!, isInitialized, isMultihypo
"""
function factorCanInitFromOtherVars(dfg::AbstractDFG,
                                    fct::Symbol,
                                    loovar::Symbol)::Tuple{Bool, Vector{Symbol}, Vector{Symbol}}
  #
  # all variables attached to this factor
  varsyms = DFG.getNeighbors(dfg, fct)

  # which element is being solved for
  sfidx = (1:length(varsyms))[varsyms .== loovar][1]
  # list of factors to use in init operation
  fctlist = []
  # list fo variables that cannot be used
  faillist = Symbol[]
  isinit = Bool[]
  for vsym in varsyms
    # check each variable one by one
    xi = DFG.getVariable(dfg, vsym)
    isi = isInitialized(xi)
    push!(isinit, isi)
    if !isi
      push!(faillist, vsym)
    end
  end

  ## determine if this factor can be used
  # priors and general n-ary cases
  canuse = length(varsyms)==1 || (length(faillist)==1 && loovar in faillist)
  ## special multihypo case (at least one hypothesis is available or initializing first hypo)
  fctnode = getFactor(dfg, fct)
  # @show canuse, isMultihypo(fctnode), isinit
  if !canuse && isMultihypo(fctnode)
    # multihypo=[1;0.5;0.5] : sfidx=1, isinit=[0,1,0] -- true
    # multihypo=[1;0.5;0.5] : sfidx=1, isinit=[0,0,1] -- true
    # multihypo=[1;0.5;0.5] : sfidx=2|3, isinit=[1,0,0] -- true
    mhp = getMultihypoDistribution(fctnode).p
    allmhp,certainidx,uncertnidx = getHypothesesVectors(mhp)
    if isLeastOneHypoAvailable(sfidx, certainidx, uncertnidx, isinit)
       # special case works
       @info "allowing init from incomplete set of previously initialized hypotheses, fct=$fct"
       canuse = true
    end
  end

  # should add the factor for use?
  if canuse
    push!(fctlist, fct)
  end

  # return if can use, the factor in an array, and the non-initialized variables attached to the factor
  return (canuse, fctlist, faillist )
end


# wow, that was quite far off -- needs testing
# function factorCanInitFromOtherVars(dfg::T,
#                                     fct::Symbol,
#                                     loovar::Symbol)::Tuple{Bool, Vector{Symbol}, Vector{Symbol}} where T <: AbstractDFG
#   #
#   # all variables attached to this factor
#   varsyms = getNeighbors(dfg, fct)
#
#   # list of factors to use in init operation
#   useinitfct = Symbol[]
#   faillist = Symbol[]
#   for vsym in varsyms
#     xi = DFG.getVariable(dfg, vsym)
#     if (isInitialized(xi) && sum(useinitfct .== fct) == 0 ) || length(varsyms) == 1
#       push!(useinitfct, fct)
#     end
#   end
#
#   return (length(useinitfct)==length(varsyms)&&length(faillist)==0,
#           useinitfct,
#           faillist   )
# end

"""
    $(SIGNATURES)

EXPERIMENTAL: initialize target variable `xi` based on connected factors in the
factor graph `fgl`.  Possibly called from `addFactor!`, or `doCliqAutoInitUp!` (?).

Notes:
- Special carve out for multihypo cases, see issue 427.

Development Notes:
> Target factor is first (singletons) or second (dim 2 pairwise) variable vertex in `xi`.
* TODO use DFG properly with local operations and DB update at end.
* TODO get faster version of `isInitialized` for database version.
* TODO: Persist this back if we want to here.
"""
function doautoinit!(dfg::T,
                     xi::DFGVariable;
                     singles::Bool=true,
                     N::Int=100,
                     logger=ConsoleLogger() )::Bool where T <: AbstractDFG
  #
  didinit = false
  # don't initialize a variable more than once
  if !isInitialized(xi)
    with_logger(logger) do
      @info "try doautoinit! of $(xi.label)"
    end
    # get factors attached to this variable xi
    vsym = xi.label
    neinodes = DFG.getNeighbors(dfg, vsym)
    # proceed if has more than one neighbor OR even if single factor
    if (singles || length(neinodes) > 1)
      # Which of the factors can be used for initialization
      useinitfct = Symbol[]
      # Consider factors connected to $vsym...
      for xifct in neinodes
        canuse, usefct, notusevars = factorCanInitFromOtherVars(dfg, xifct, vsym)
        if canuse
          union!(useinitfct, usefct)
        end
      end
      with_logger(logger) do
        @info "init with useinitfct $useinitfct"
      end
      # println("Consider all singleton (unary) factors to $vsym...")
      # calculate the predicted belief over $vsym
      if length(useinitfct) > 0
        with_logger(logger) do
          @info "do init of $vsym"
        end
        pts,inferdim = predictbelief(dfg, vsym, useinitfct, logger=logger)
        setValKDE!(xi, pts, true, inferdim)
        # Update the estimates (longer DFG function used so cloud is also updated)
        setVariablePosteriorEstimates!(dfg, xi.label)
        # Update the data in the event that it's not local
        updateVariableSolverData!(dfg, xi, :default, true)    # TODO perhaps usecopy=false
        # deepcopy graphinit value, see IIF #612
        updateVariableSolverData!(dfg, xi.label, getSolverData(xi, :default), :graphinit, true, Symbol[], false)
        didinit = true
      end
    end
  end
  return didinit
end

function doautoinit!(dfg::T,
                     Xi::Vector{<:DFGVariable};
                     singles::Bool=true,
                     N::Int=100,
                     logger=ConsoleLogger() )::Bool where T <: AbstractDFG
  #
  #
  # Mighty inefficient function, since we only need very select fields nearby from a few neighboring nodes
  # do double depth search for variable nodes

  didinit = true

  # loop over all requested variables that must be initialized
  for xi in Xi
    didinit &= doautoinit!(dfg, xi, singles=singles, N=N, logger=logger)
  end
  return didinit
end

function doautoinit!(dfg::T,
                     xsyms::Vector{Symbol};
                     singles::Bool=true,
                     N::Int=100,
                     logger=SimpleLogger(logger)  )::Bool where T <: AbstractDFG
  #
  verts = getVariable.(dfg, xsyms)
  return doautoinit!(dfg, verts, singles=singles, N=N, logger=logger)
end
function doautoinit!(dfg::T,
                     xsym::Symbol;
                     singles::Bool=true,
                     N::Int=100,
                     logger=ConsoleLogger()  )::Bool where T <: AbstractDFG
  #
  return doautoinit!(dfg, [getVariable(dfg, xsym);], singles=singles, N=N, logger=logger)
end

"""
    $(SIGNATURES)

Workaround function when first-version (factor graph based) auto initialization fails.  Usually occurs when using factors that have high connectivity to multiple variables.
"""
function initManual!(dfg::AbstractDFG, vert::DFGVariable, pX::BallTreeDensity)::Nothing
  setValKDE!(vert, pX, true)
  return nothing
end
function initManual!(dfg::AbstractDFG, sym::Symbol, pX::BallTreeDensity)::Nothing
  vert = getVariable(dfg, sym)
  initManual!(dfg, vert, pX)
  return nothing
end
function initManual!(dfg::AbstractDFG, sym::Symbol, usefcts::Vector{Symbol})::Nothing
  @info "initManual! $sym"
  pts = predictbelief(dfg, sym, usefcts)
  vert = getVariable(dfg, sym)
  Xpre = AMP.manikde!(pts, getSofttype(vert) |> getManifolds )
  setValKDE!(vert, Xpre, true)
  return nothing
end


function initManual!(dfg::AbstractDFG, sym::Symbol, pts::Array{Float64,2})
  var = getVariable(dfg, sym)
  pp = manikde!(pts, getManifolds(var))
  initManual!(dfg,sym,pp)
end

const initVariableManual! = initManual!

"""
    $SIGNATURES

Set solveKey values of `dest::AbstractDFG` according to `initKey::Symbol=:graphinit` values.

Notes
- Some flexibility for using two DFGs and different key values, see Examples and code for details.
- Can also be specific with `varList::Vector{Symbol}`.
- Returns `dest` graph.
- Uses the supersolve mechanism.

Examples
```julia
resetInitialValues!(fg)
resetInitialValues!(fg1,fg2)  # into 1 from 2
resetInitialValues!(fg1,fg1,:myotherinit)  # use different init value into solveKey :default
resetInitialValues!(fg1,fg1,:graphinit, :mysolver) # not into solveKey=:default but :mysolver
resetInitialValues!(fg1, varList=[:x1;:l3])  # Specific variables only

# Into `fgNew` object, leaving `fg` untouched
fgNew = deepcopy(fg)
resetInitialValues!(fgNew,fg)
```

Related

initManual!, graphinit (keyword)
"""
function resetInitialValues!(dest::AbstractDFG,
                             src::AbstractDFG=dest,
                             initKey::Symbol=:graphinit,
                             solveKey::Symbol=:default;
                             varList::AbstractVector{Symbol}=ls(dest))
  #
  for vs in varList
    vnd = getSolverData(getVariable(src, vs), initKey)
    # guess we definitely want to use copy to preserve the initKey memory
    updateVariableSolverData!(dest,vs,vnd,solveKey,true)
  end
  return dest
end
const resetInitValues! = resetInitialValues!

"""
    $SIGNATURES

Ensure that no variables set as `solvable=1` are floating free without any connected `solvable=1` factors.  If any found, then set those 'free' variable's `solvable=solvableFallback` (default `0`).

Related

ensureAllInitialized!
"""
function ensureSolvable!(dfg::AbstractDFG; solvableTarget::Int=1, solvableFallback::Int=0)
  # workaround in case isolated variables occur
  solvVars = ls(dfg, solvable=solvableTarget)
  varHasFact = (x->length(ls(dfg,x, solvable=solvableTarget))==0).(solvVars)
  blankVars = solvVars[findall(varHasFact)]
  if 0 < length(blankVars)
    @warn("solveTree! dissallows solvable variables without any connected solvable factors -- forcing solvable=0 on $(blankVars)")
    (x->setSolvable!(dfg, x, solvableFallback)).(blankVars)
  end
  return blankVars
end

"""
    $SIGNATURES

Perform `graphinit` over all variables with `solvable=1` (default).

Related

ensureSolvable!, (EXPERIMENTAL 'treeinit')
"""
function ensureAllInitialized!(dfg::T; solvable::Int=1) where T <: AbstractDFG
  # allvarnodes = getVariables(dfg)
  syms = intersect(getAddHistory(dfg), ls(dfg, solvable=solvable) )
  # syms = ls(dfg, solvable=solvable) # |> sortDFG
  repeatCount = 0
  repeatFlag = true
  while repeatFlag
    repeatFlag = false
    repeatCount += 1
    if 10 < repeatCount
      @info "not able to initialize all variables via the factor graph, abort autoinit."
      break;
    end
    for sym in syms
      var = getVariable(dfg, sym)
      if !isInitialized(var)
        @info "$(var.label) is not initialized, and will do so now..."
        doautoinit!(dfg, [var;], singles=true)
        !isInitialized(var) ? (repeatFlag = true) : nothing
      end
    end
  end
  nothing
end

function assembleFactorName(dfg::AbstractDFG,
                            Xi::Vector{<:DFGVariable};
                            maxparallel::Union{Nothing,Int}=nothing)
  #
  if maxparallel !== nothing
    @warn "keyword maxparallel has been deprecated, use getSolverParams(fg).maxincidence=$maxparallel instead."
    getSolverParams(dfg).maxincidence = maxparallel
  end

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
"""
function addFactor!(dfg::AbstractDFG,
                    Xi::Vector{<:DFGVariable},
                    usrfnc::R;
                    multihypo::Vector{Float64}=Float64[],
                    nullhypo::Float64=0.0,
                    solvable::Int=1,
                    tags::Vector{Symbol}=Symbol[],
                    timestamp::Union{DateTime,ZonedDateTime}=now(localzone()),
                    graphinit::Bool=getSolverParams(dfg).graphinit,
                    threadmodel=SingleThreaded,
                    maxparallel::Union{Int,Nothing}=nothing  ) where
                      {R <: FunctorInferenceType}
  #
  # depcrecation
  if maxparallel !== nothing
    @warn "maxparallel keyword is deprecated, use getSolverParams(fg).maxincidence instead."
    getSolverParams(dfg).maxincidence = maxparallel
  end

  varOrderLabels = [v.label for v=Xi]
  namestring = assembleFactorName(dfg, Xi)
  solverData = getDefaultFactorData(dfg, Xi, deepcopy(usrfnc), multihypo=multihypo, nullhypo=nullhypo, threadmodel=threadmodel)
  newFactor = DFGFactor(Symbol(namestring),
                        varOrderLabels,
                        solverData;
                        tags=Set(union(tags, [:FACTOR])),
                        solvable=solvable,
                        timestamp=timestamp)
  #

  success = DFG.addFactor!(dfg, newFactor)

  # TODO: change this operation to update a conditioning variable
  graphinit && doautoinit!(dfg, Xi, singles=false)

  return newFactor
end

function addFactor!(dfg::AbstractDFG,
                    xisyms::Vector{Symbol},
                    usrfnc::FunctorInferenceType;
                    multihypo::Vector{<:Real}=Float64[],
                    nullhypo::Float64=0.0,
                    solvable::Int=1,
                    timestamp::Union{DateTime,ZonedDateTime}=now(localzone()),
                    tags::Vector{Symbol}=Symbol[],
                    graphinit::Bool=getSolverParams(dfg).graphinit,
                    threadmodel=SingleThreaded,
                    maxparallel::Union{Nothing,Int}=nothing  )
  #
  # depcrecation
  if maxparallel !== nothing
    @warn "maxparallel keyword is deprecated, use getSolverParams(fg).maxincidence instead."
    getSolverParams(dfg).maxincidence = maxparallel
  end

  verts = map(vid -> DFG.getVariable(dfg, vid), xisyms)
  addFactor!(dfg, verts, usrfnc, multihypo=multihypo, nullhypo=nullhypo, solvable=solvable, tags=tags, graphinit=graphinit, threadmodel=threadmodel, timestamp=timestamp )
end




function prtslperr(s)
  println(s)
  sleep(0.1)
  error(s)
end

"""
    $SIGNATURES

Determine the variable ordering used to construct both the Bayes Net and Bayes/Junction/Elimination tree.

Notes
- Heuristic method -- equivalent to QR or Cholesky.
- Are using Blas `QR` function to extract variable ordering.
- **NOT USING SUITE SPARSE** -- which would requires commercial license.
- For now `A::Array{<:Number,2}` as a dense matrix.
- Columns of `A` are system variables, rows are factors (without differentiating between partial or full factor).
- default is to use `solvable=1` and ignore factors and variables that might be used for dead reckoning or similar.

Future
- TODO: `A` should be sparse data structure (when we exceed 10'000 var dims)
- TODO: Incidence matrix is rectagular and adjacency is the square.
"""
function getEliminationOrder(dfg::G;
                             ordering::Symbol=:qr,
                             solvable::Int=1,
                             constraints::Vector{Symbol}=Symbol[]) where G <: AbstractDFG
  #
  @assert 0 == length(constraints) || ordering == :ccolamd "Must use ordering=:ccolamd when trying to use constraints"
  # Get the sparse adjacency matrix, variable, and factor labels
  adjMat, permuteds, permutedsf = DFG.getBiadjacencyMatrix(dfg, solvable=solvable)
  # adjMat, permuteds, permutedsf = DFG.getAdjacencyMatrixSparse(dfg, solvable=solvable)

  # Create dense adjacency matrix

  p = Int[]
  if ordering==:chol
    # hack for dense matrix....
    A = Array(adjMat)
    p = cholfact(A'A,:U,Val(true))[:p] #,pivot=true
    @warn "check cholesky ordering is not reversed -- basically how much fill in (separator size) are you seeing???  Long skinny chains in tree is bad."
  elseif ordering==:qr
    # hack for dense matrix....
    A = Array(adjMat)
    # this is the default
    q,r,p = qr(A, Val(true))
    p .= p |> reverse
  elseif ordering==:ccolamd
    cons = zeros(SuiteSparse_long, length(adjMat.colptr) - 1)
    cons[findall(x->x in constraints, permuteds)] .= 1
    p = Ccolamd.ccolamd(adjMat, cons)
    @warn "Ccolamd is experimental in IIF at this point in time."
  else
    prtslperr("getEliminationOrder -- cannot do the requested ordering $(ordering)")
  end

  # Return the variable ordering that we should use for the Bayes map
  # reverse order checked in #475 and #499
  return permuteds[p]
end


# lets create all the vertices first and then deal with the elimination variables thereafter
function addBayesNetVerts!(dfg::AbstractDFG,
                           elimOrder::Array{Symbol,1} )
  #
  for pId in elimOrder
    vert = DFG.getVariable(dfg, pId)
    if getSolverData(vert).BayesNetVertID == nothing || getSolverData(vert).BayesNetVertID == :_null # Special serialization case of nothing
      @debug "[AddBayesNetVerts] Assigning $pId.data.BayesNetVertID = $pId"
      getSolverData(vert).BayesNetVertID = pId
    else
      @warn "addBayesNetVerts -- Something is wrong, variable '$pId' should not have an existing Bayes net reference to '$(getSolverData(vert).BayesNetVertID)'"
    end
  end
end

function addConditional!(dfg::AbstractDFG,
                         vertId::Symbol,
                         Si::Vector{Symbol} )
  #
  bnv = DFG.getVariable(dfg, vertId)
  bnvd = getSolverData(bnv)
  bnvd.separator = Si
  for s in Si
    push!(bnvd.BayesNetOutVertIDs, s)
  end
  return nothing
end

function addChainRuleMarginal!(dfg::AbstractDFG,
                               Si::Vector{Symbol};
                               maxparallel::Union{Nothing,Int}=nothing )
  #
  if maxparallel !== nothing
    @warn "keyword maxparallel has been deprecated, use getSolverParams(fg).maxincidence=$maxparallel instead."
    getSolverParams(dfg).maxincidence = maxparallel
  end

  lbls = String[]
  genmarg = GenericMarginal()
  Xi = map(v -> DFG.getVariable(dfg, v), Si)
  # @info "adding marginal to"
  # for x in Xi
  #   @info "x.index=",x.index
  # end
  addFactor!( dfg, Xi, genmarg, graphinit=false )
  nothing
end

function rmVarFromMarg(dfg::AbstractDFG,
                       fromvert::DFGVariable,
                       gm::Vector{DFGFactor};
                       maxparallel::Union{Nothing,Int}=nothing  )
  #
  if maxparallel !== nothing
    @warn "keyword maxparallel has been deprecated, use getSolverParams(fg).maxincidence=$maxparallel instead."
    getSolverParams(dfg).maxincidence = maxparallel
  end
  @info " - Removing $(fromvert.label)"
  for m in gm
    @info "Looking at $(m.label)"
    for n in DFG.getNeighbors(dfg, m) #x1, x2
      if n == fromvert.label # n.label ==? x1
        @info "   - Breaking link $(m.label)->$(fromvert.label)..."
        @info "     - Original links: $(DFG.ls(dfg, m))"
        remvars = setdiff(DFG.ls(dfg, m), [fromvert.label])
        @info "     - New links: $remvars"

        DFG.deleteFactor!(dfg, m) # Remove it
        if length(remvars) > 0
          @info "$(m.label) still has links to other variables, readding it back..."
          addFactor!(dfg, remvars, getSolverData(m).fnc.usrfnc!, graphinit=false )
        else
          @info "$(m.label) doesn't have any other links, not adding it back..."
        end
      end
    end
    # Added back in chain rule.
    if DFG.exists(dfg, m) && length(DFG.getNeighbors(dfg, m)) <= 1
      @warn "removing vertex id=$(m.label)"
      DFG.deleteFactor!(dfg, m)
    end
  end
  return nothing
end

function buildBayesNet!(dfg::AbstractDFG,
                        elimorder::Vector{Symbol};
                        maxparallel::Union{Nothing,Int}=nothing,
                        solvable::Int=1 )
  #
  # addBayesNetVerts!(dfg, elimorder)
  for v in elimorder
    @info ""
    @info "Eliminating $(v)"
    @info "==============="
    @info ""
    # which variable are we eliminating

    # all factors adjacent to this variable
    fi = Symbol[]
    Si = Symbol[]
    gm = DFGFactor[]

    vert = DFG.getVariable(dfg, v)
    for fctId in DFG.getNeighbors(dfg, vert, solvable=solvable)
      fct = DFG.getFactor(dfg, fctId)
      if (getSolverData(fct).eliminated != true)
        push!(fi, fctId)
        for sepNode in DFG.getNeighbors(dfg, fct, solvable=solvable)
          # TODO -- validate !(sepNode.index in Si) vs. older !(sepNode in Si)
          if sepNode != v && !(sepNode in Si) # Symbol comparison!
            push!(Si,sepNode)
          end
        end
        getSolverData(fct).eliminated = true
      end

      if typeof(getSolverData(fct).fnc) == CommonConvWrapper{GenericMarginal}
        push!(gm, fct)
      end
    end

    if v != elimorder[end]
      addConditional!(dfg, v, Si)
      # not yet inserting the new prior p(Si) back into the factor graph
    end

    # mark variable
    getSolverData(vert).eliminated = true

    # TODO -- remove links from current vertex to any marginals
    rmVarFromMarg(dfg, vert, gm)

    #add marginal on remaining variables... ? f(xyz) = f(x | yz) f(yz)
    # new function between all Si (round the outside, right the outside)
    length(Si) > 0 && addChainRuleMarginal!(dfg, Si)

  end
  return nothing
end



"""
    $(SIGNATURES)

Get KernelDensityEstimate kde estimate stored in variable node.
"""
getBelief(vnd::VariableNodeData) = AMP.manikde!(getVal(vnd), getBW(vnd)[:,1], getSofttype(vnd) |> getManifolds)
getBelief(v::DFGVariable, solvekey::Symbol=:default) = getKDE(getSolverData(v, solvekey))
getBelief(dfg::AbstractDFG, lbl::Symbol, solvekey::Symbol=:default) = getKDE(getVariable(dfg, lbl), solvekey)

const getKDE = getBelief


#
