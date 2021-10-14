
##==============================================================================
## LEGACY SUPPORT FOR ZMQ IN CAESAR
##==============================================================================

export listSolvekeys

export _evalType

# not sure if and where this is still being used
function _evalType(pt::String)::Type
  @error "_evalType has been deprecated, use DFG serialization methods instead."
  try
    getfield(Main, Symbol(pt))
  catch ex
    io = IOBuffer()
    showerror(io, ex, catch_backtrace())
    err = String(take!(io))
    error("_evalType: Unable to locate factor/distribution type '$pt' in main context (e.g. do a using on all your factor libraries). Please check that this factor type is loaded into main. Stack trace = $err")
  end
end


##==============================================================================
## Deprecate code below before v0.27
##==============================================================================


getOutNeighbors(w...;kw...) = error("Obsolete, use DFG.getNeighbors instead.")

# # TODO -- there should be a better way, without retrieving full vertex
# # TODO -- Deprecated for DFG -- must update
# function getOutNeighbors(dfg::T, v::V; needdata::Bool=false, ready::Union{Nothing, Int}=nothing, backendset::Union{Nothing, Int}=nothing)::Vector{Symbol} where {T <: AbstractDFG, V <: DFGNode}
#   @warn "TODO: needdata is currently ignored. Symbols are returned."
#   nodes = getNeighbors(dfg, v, ready=ready, backendset=backendset)
#   return nodes
# end
# function getOutNeighbors(dfg::T, vertSym::Symbol; needdata::Bool=false, ready::Int=1, backendset::Int=1 )::Vector{Symbol} where {T <: AbstractDFG, V <: DFGNode}
#   @warn "TODO: needdata is currently ignored. Symbols are returned."
#   nodes = getNeighbors(dfg, vertSym, ready=ready, backendset=backendset)
#   return nodes
# end

# SolverParams(;dimID::Int=0,
#               registeredModuleFunctions=nothing,
#               reference=nothing,
#               stateless::Bool=false,
#               qfl::Int=99999999999,
#               isfixedlag::Bool=false,
#               limitfixeddown::Bool=false,
#               incremental::Bool=true,
#               useMsgLikelihoods::Bool=false,
#               upsolve::Bool=true,
#               downsolve::Bool=true,
#               drawtree::Bool=false,
#               drawCSMIters::Bool=true,
#               showtree::Bool=false,
#               drawtreerate::Float64=0.5,
#               dbg::Bool=false,
#               async::Bool=false,
#               limititers::Int=500,
#               N::Int=100,
#               multiproc::Bool=1 < nprocs(),
#               logpath::String="/tmp/caesar/$(now())",
#               graphinit::Bool=true,
#               treeinit::Bool=false,
#               limittreeinit_iters::Int=10,
#               algorithms::Vector{Symbol}=[:default],
#               spreadNH::Real=3.0,
#               inflation::Real=5.0,
#               inflateCycles::Int=3,
#               gibbsIters::Int=3,
#               maxincidence::Int=500,
#               alwaysFreshMeasurements::Bool=true,
#               attemptGradients::Bool=true,
#               devParams::Dict{Symbol,String}=Dict{Symbol,String}()
#             ) = begin useMsgLikelihoods==true && @warn "useMsgLikelihoods is under development, use with care, see #1010"
#                 SolverParams( dimID,
#                               registeredModuleFunctions,
#                               reference,
#                               stateless,
#                               qfl,
#                               isfixedlag,
#                               limitfixeddown,
#                               incremental,
#                               useMsgLikelihoods,
#                               upsolve,
#                               downsolve,
#                               drawtree,
#                               drawCSMIters,
#                               showtree,
#                               drawtreerate,
#                               dbg,
#                               async,
#                               limititers,
#                               N,
#                               multiproc,
#                               logpath,
#                               graphinit,
#                               treeinit,
#                               limittreeinit_iters,
#                               algorithms,
#                               spreadNH,
#                               inflation,
#                               inflateCycles,
#                               gibbsIters,
#                               maxincidence,
#                               alwaysFreshMeasurements,
#                               attemptGradients,
#                               devParams )
#             end
#

@deprecate setVariableInferDim!(w...;kw...) setIPC!(w...;kw...)

# moved to IncrementalInference/attic
# include("CSMOccuranceUtils.jl")

@deprecate prepgenericconvolution(w...;kw...) _prepCCW(w...;kw...)

function getSample(cf::CalcFactor, N::Int) 
  # Base.depwarn("`getSample(cf, N)` is deprecated, use `getSample(cf)`", :getSample)
  error("`getSample(cf, N)` is deprecated, use `getSample(cf)`")
end

# """
#     $SIGNATURES

# Build an approximate density `[Y|X,DX,.]=[X|Y,DX][DX|.]` as proposed by the conditional convolution.

# Notes
# - Assume both are on circular manifold, `manikde!(pts, (:Circular,))`
# """
# function approxConvCircular(pX::ManifoldKernelDensity, 
#                             pDX::ManifoldKernelDensity; N::Int=100)
#   #

#   # building basic factor graph
#   tfg = initfg()
#   addVariable!(tfg, :s1, Sphere1)
#   addVariable!(tfg, :s2, Sphere1)
#   addFactor!(tfg, [:s1;:s2], Sphere1Sphere1(pDX), graphinit=false)
#   initManual!(tfg,:s1, pX)

#   # solve for outgoing proposal value
#   approxConv(tfg,:s1s2f1,:s2)
# end

# function approxConvCircular(pX::ManifoldKernelDensity, 
#                             pDX::SamplableBelief; N::Int=100)
#   #
#   pts = reshape(rand(pDX, N), 1, :)
#   pC = manikde!(pts, Sphere1)
#   approxConvCircular(pX, pC)
# end


# function approxConvCircular(pX::SamplableBelief, 
#                             pDX::ManifoldKernelDensity; N::Int=100)
#   #
#     pts = reshape(rand(pX, N), 1, :)
#   pC = manikde!(pts, Sphere1)
#   approxConvCircular(pC, pDX)
# end



@deprecate printGraphSummary(dfg::AbstractDFG, logger=ConsoleLogger()) show(logger.stream, MIME("text/plain"), dfg)
@deprecate printSummary(dfg::AbstractDFG, logger=ConsoleLogger()) show(logger.stream, MIME("text/plain"), dfg)


# """
#     $SIGNATURES

# Print basic summary of graph to `logger=ConsoleLogger()`.
# """
# function printGraphSummary(dfg::G, logger=ConsoleLogger())::Nothing where {G <: AbstractDFG}
#     vars = ls(dfg)
#     fcts = lsf(dfg)

#     prio = lsfPriors(dfg)

#     isinit = map(x->isInitialized(dfg,x), vars)
#     infdim = map(x->getVariableInferredDim(dfg, x), vars)
#     numedges = map(v->length(ls(dfg, v)), vars)
#     numfed = map(fc->length(ls(dfg, fc)), fcts)
#     vardims = map(v->getDimension(getVariable(dfg, v)), vars)
#     fctdims = map(v->getDimension(getFactor(dfg, v)), fcts)
#     priodims = map(v->getDimension(getFactor(dfg, v)), prio)

#     with_logger(logger) do
#       @info "Distributed Factor Graph summary:"
#       @info "  num variables:    $(length(vars))"
#       @info "  num factors:      $(length(fcts)), w/ $(length(prio)) priors"
#       @info "  var initialized:  $(sum(isinit))"
#       @info ""
#       @info "  var num edges: min. $(minimum(numedges)) | mean $(round(Statistics.mean(numedges),digits=2)) | 90% $(round(quantile(numedges,0.9),digits=2)) | max. $(maximum(numedges))"
#       @info "  fct num edges: min. $(minimum(numfed)) | mean $(round(Statistics.mean(numfed),digits=2)) | 90% $(round(quantile(numfed,0.9),digits=2)) | max. $(maximum(numfed))"
#       @info "  Variable dims: min. $(minimum(vardims)) | mean $(round(Statistics.mean(vardims),digits=2)) | 90% $(round(quantile(vardims,0.9),digits=2)) | max. $(maximum(vardims))"
#       @info "  Factor dims:   min. $(minimum(fctdims)) | mean $(round(Statistics.mean(fctdims),digits=2)) | 90% $(round(quantile(fctdims,0.9),digits=2)) | max. $(maximum(fctdims))"
#       @info "  Prior dimens:  min. $(minimum(priodims)) | mean $(round(Statistics.mean(priodims),digits=2)) | 90% $(round(quantile(priodims,0.9),digits=2)) | max. $(maximum(priodims))"
#       @info "  var infr'dims: min. $(minimum(infdim)) | mean $(round(Statistics.mean(infdim),digits=2)) | 90% $(round(quantile(infdim,0.9),digits=2)) | max. $(maximum(infdim))"
#     end
#     nothing
# end

# """
#     $SIGNATURES

# Print basic summary of graph to `logger=ConsoleLogger()`.
# """
# function printSummary(dfg::G, logger=ConsoleLogger()) where G <: AbstractDFG
#     printGraphSummary(dfg, logger)
# end



# """
#     $SIGNATURES

# Calculate both measured and predicted relative variable values, starting with `from` at zeros up to `to::Symbol`.

# Notes
# - assume single variable separators only.

# DevNotes
# - TODO better consolidate with [`approxConvBelief`](@ref) which can now also work with factor chains.
# """
# function accumulateFactorChain( dfg::AbstractDFG,
#                                 from::Symbol,
#                                 to::Symbol,
#                                 fsyms::Vector{Symbol}=findFactorsBetweenNaive(dfg, from, to);
#                                 initval=zeros(size(getVal(dfg, from))))

#   # get associated variables
#   svars = union(ls.(dfg, fsyms)...)

#   # use subgraph copys to do calculations
#   tfg_meas = buildSubgraph(dfg, [svars;fsyms])
#   tfg_pred = buildSubgraph(dfg, [svars;fsyms])

#   # drive variable values manually to ensure no additional stochastics are introduced.
#   nextvar = from
#   initManual!(tfg_meas, nextvar, initval)
#   initManual!(tfg_pred, nextvar, initval)

#   # nextfct = fsyms[1] # for debugging
#   for nextfct in fsyms
#     nextvars = setdiff(ls(tfg_meas,nextfct),[nextvar])
#     @assert length(nextvars) == 1 "accumulateFactorChain requires each factor pair to separated by a single variable"
#     nextvar = nextvars[1]
#     meas, pred = approxDeconv(dfg, nextfct) # solveFactorMeasurements
#     pts_meas = approxConv(tfg_meas, nextfct, nextvar, (meas,ones(Int,100),collect(1:100)))
#     pts_pred = approxConv(tfg_pred, nextfct, nextvar, (pred,ones(Int,100),collect(1:100)))
#     initManual!(tfg_meas, nextvar, pts_meas)
#     initManual!(tfg_pred, nextvar, pts_pred)
#   end
#   return getVal(tfg_meas,nextvar), getVal(tfg_pred,nextvar)
# end


# # TODO should this be consolidated with regular approxConv?
# # TODO, perhaps pass Xi::Vector{DFGVariable} instead?
# function approxConvBinary(arr::Vector{Vector{Float64}},
#                           meas::AbstractFactor,
#                           outdims::Int,
#                           fmd::FactorMetadata,
#                           measurement::Tuple=(Vector{Vector{Float64}}(),);
#                           varidx::Int=2,
#                           N::Int=length(arr),
#                           vnds=DFGVariable[],
#                           _slack=nothing )
#   #
#   # N = N == 0 ? size(arr,2) : N
#   pts = [zeros(outdims) for _ in 1:N];
#   ptsArr = Vector{Vector{Vector{Float64}}}()
#   push!(ptsArr,arr)
#   push!(ptsArr,pts)

#   fmd.arrRef = ptsArr

#   # TODO consolidate with ccwl??
#   # FIXME do not divert Mixture for sampling
#   # cf = _buildCalcFactorMixture(ccwl, fmd, 1, ccwl.measurement, ARR) # TODO perhaps 0 is safer
#   # FIXME 0, 0, ()
#   cf = CalcFactor( meas, fmd, 0, 0, (), ptsArr)

#   measurement = length(measurement[1]) == 0 ? sampleFactor(cf, N) : measurement
#   # measurement = size(measurement[1],2) == 0 ? sampleFactor(meas, N, fmd, vnds) : measurement

#   zDim = length(measurement[1][1])
#   ccw = CommonConvWrapper(meas, ptsArr[varidx], zDim, ptsArr, fmd, varidx=varidx, measurement=measurement)  # N=> size(measurement[1],2)

#   for n in 1:N
#     ccw.cpt[Threads.threadid()].particleidx = n
#     _solveCCWNumeric!( ccw, _slack=_slack )
#   end
#   return pts
# end

@deprecate getParametricMeasurement(w...;kw...) getMeasurementParametric(w...;kw...)

# function prtslperr(s)
#   println(s)
#   sleep(0.1)
#   error(s)
# end

@deprecate getKDE(w...;kw...) getBelief(w...;kw...)

@deprecate findRelatedFromPotential(w...;kw...) (calcProposalBelief(w...;kw...),)

# function generateNullhypoEntropy( val::AbstractMatrix{<:Real},
#                                   maxlen::Int,
#                                   spreadfactor::Real=10  )
#   #
#   # covD = sqrt.(vec(Statistics.var(val,dims=2))) .+ 1e-3
#   # cVar = diagm((spreadfactor*covD).^2)
#   len = size(val, 1)
#   cVar = diagm((spreadfactor*ones(len)).^2)
#   mu = zeros( len )
#   MvNormal( mu, cVar )
# end

# @deprecate calcFactorResidualTemporary( fct::AbstractRelative, measurement::Tuple, T_param_args...; kw...) calcFactorResidualTemporary(fct, measurement, tuple(((x->x[1]).(T_param_args))...), tuple(((x->x[2]).(T_param_args))...); kw...)

# @deprecate _buildGraphByFactorAndTypes!(fct::AbstractFactor, TypeParams_vec...; kw...) _buildGraphByFactorAndTypes!(fct, (x->x[1]).(TypeParams_vec), (x->x[2]).(TypeParams_vec) ; kw...)

# @deprecate _evalFactorTemporary!( fct::AbstractFactor, sfidx::Int, meas_single::Tuple, TypeParams_args...; kw...) _evalFactorTemporary!( fct,sfidx,meas_single, (x->x[1]).(TypeParams_args), (x->x[2]).(TypeParams_args); kw...)

# MOVED TO AMP v0.4.6
# # NOTE SWITCHED TO ON-MANIFOLD VERSION, but this used to give deviation
# calcVariableCovarianceBasic(M::AbstractManifold, ptsArr::Vector{Vector{Float64}})
#   # cannot calculate the stdev from uninitialized state
#   # FIXME assume point type is only Vector{Float} at this time
#   @cast arr[i,j] := ptsArr[j][i]
#   msst = Statistics.std(arr, dims=2)
#   # FIXME use adaptive scale, see #802
#   msst_ = 0 < sum(1e-10 .< msst) ? maximum(msst) : 1.0
#   return msst_
# end

@deprecate testFactorResidualBinary(fct, T1::InstanceType{InferenceVariable}, T2::InstanceType{InferenceVariable}, param1, param2, meas) calcFactorResidualTemporary(fct, (T1,T2), meas, (param1,param2))

# """
#     $SIGNATURES
#
# Return the directly achievable dimensionality of solve for each variable in a clique.
#
# Related
#
# getFactorSolvableDim
# """
# function getCliqVarPossibleDim(dfg::G,
#                                     cliq::TreeClique,
#                                     saturate::Bool=true,
#                                     fraction::Bool=true  )::Dict{Symbol, Float64}
#   #
#   # variables and factors associated with this clique
#   vars = getCliqAllVarIds(cliq)
#   # fcts = getCliqAllFactIds(cliq)
#   # rows = length(fcts)
#   cols = length(vars)
#
#   # for output result
#   dict = Dict{Symbol,Float64}()
#
#   for j in 1:cols
#     dict[vars[j]] += getVariablePossibleDim(dfg, vars[j])
#   end
#
# end


# getManifold(::InstanceType{LinearRelative{N}}) where {N} = Euclidean(N)
# getManifolds(::T) where {T <: LinearRelative} = convert(Tuple, getManifold(T))
# getManifolds(::Type{<:T}) where {T <: LinearRelative} = convert(Tuple, getManifold(T))
# getManifolds(fctType::Type{LinearRelative}) = getManifolds(getDomain(fctType))

# # FIXME, why is Manifolds depdendent on the solveKey?? Should just be at DFGVariable level?

# getManifolds(vd::VariableNodeData) = getVariableType(vd) |> getManifolds
# getManifolds(::DFGVariable{T}) where T <: InferenceVariable = getManifolds(T)
# getManifolds(dfg::AbstractDFG, sym::Symbol) = getManifolds(getVariable(dfg, sym))
#
# getManifolds(vartype::InferenceVariable) = vartype.manifolds
# getManifolds(vartype::Type{<: InferenceVariable}) = getManifolds(vartype())


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

# function manikde!(ptsArr::Vector{Vector{Float64}}, manis::Tuple)
#   arr = Matrix{Float64}(undef, length(ptsArr[1]), length(ptsArr))
#   @cast arr[i,j] = ptsArr[j][i]
#   manikde!(arr, manis)
# end


##==============================================================================
## Deprecate code below before v0.26
##==============================================================================

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

@deprecate calcPerturbationFromVariable( fgc::FactorGradientsCached!, fromVar::Int, infoPerCoord::AbstractVector;tol::Real=0.02*fgc._h ) calcPerturbationFromVariable(fgc, [fromVar => infoPerCoord;], tol=tol )

function getVal(vA::Vector{<:DFGVariable}, solveKey::Symbol=:default)
  @error "getVal(::Vector{DFGVariable}) is obsolete, use getVal.(DFGVariable) instead."
  # len = length(vA)
  # vals = Array{Array{Float64,2},1}()
  # cols = Array{Int,1}()
  # push!(cols,0)
  # rows = Array{Int,1}()
  # for v in vA
  #     push!(vals, getVal(v, solveKey=solveKey))
  #     c = size(vals[end],2)
  #     r = size(vals[end],1)
  #     push!(cols, floor(Int,c))
  #     push!(rows, floor(Int,r))
  # end
  # cols = cumsum(cols)
  # sc = cols[end]
  # rw = floor(Int,rows[1])
  # val = Array{Float64,2}(undef,rw, sc)
  # for i in 1:(len-1)
  #     val[:,(cols[i]+1):cols[i+1]] = vals[i]
  # end
  # val[:,(cols[len]+1):cols[len+1]] = vals[len] # and the last one
  # return val
end

# """
#     $SIGNATURES

# Get graph node (variable or factor) dimension.
# """
# DFG.getDimension(vartype::InferenceVariable) = vartype.dims #TODO Deprecate
# DFG.getDimension(vartype::Type{<:InferenceVariable}) = getDimension(vartype())


@deprecate ensureAllInitialized!(w...;kw...) initAll!(w...;kw...)

@deprecate getFactorMean(w...) IIF.getMeasurementParametric(w...)[1]

# """
#     $SIGNATURES

# Recover the mean (Gaussian) or estimate stochastic mean (non-Gaussian) value stored in a factor measurement.

# Related

# accumulateFactorMeans, solveBinaryFactorParameteric
# """
# function getFactorMean(fct::FunctorInferenceType)
#   fctt = typeof(fct)
#   error("no getFactorMean defined for $(fctt.name), has fields $(fieldnames(fctt))")
# end

# getFactorMean(fct::Normal) = [fct.μ]
# getFactorMean(fct::MvNormal) = fct.μ
# getFactorMean(fct::Union{<:BallTreeDensity,<:ManifoldKernelDensity}) = getKDEMean(fct)
# getFactorMean(fct::AliasingScalarSampler) = Statistics.mean(rand(fct,1000))

# getFactorMean(fct::DFGFactor) = getFactorMean(getFactorType(fct))

# getFactorMean(dfg::AbstractDFG, fctsym::Symbol) = getFactorMean(getFactor(dfg, fctsym))


# AMP.getManifolds(::T) where {T <: InferenceVariable} = getManifolds(getManifold(T))
# AMP.getManifolds(::Type{T}) where {T <: InferenceVariable} = getManifolds(getManifold(T))


# getManifolds(getManifold(T))
getManifolds(::InstanceType{T}) where {T <: Union{InferenceVariable, AbstractFactor}} = error("getManifolds is obsolete, use getManifold(...)::MB.AbstactManifold instead.")





#
