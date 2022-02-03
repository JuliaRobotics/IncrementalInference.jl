
##==============================================================================
## LEGACY, towards Sidecar
##==============================================================================

"""
Converter: Prior -> PackedPrior::Dict{String, Any}

FIXME see DFG #590 for consolidation with Serialization and Marshaling
"""
function convert(::Type{Dict{String, Any}}, prior::IncrementalInference.Prior)
    @error("Obsolete, use pack/unpack converters instead")
    z = convert(Type{Dict{String, Any}}, prior.Z)
    return Packed_Factor([z], "Prior")
end

"""
Converter: PackedPrior::Dict{String, Any} -> Prior

FIXME see DFG #590 for consolidation on Serialization and Marshaling
"""
function convert(::Type{<:Prior}, prior::Dict{String, Any})
    @error("Obsolete, use pack/unpack converters instead")
    # Genericize to any packed type next.
    z = prior["measurement"][1]
    z = convert(DFG.getTypeFromSerializationModule(z["distType"]), z)
    return Prior(z)
end


##==============================================================================
## Deprecate code below before v0.29
##==============================================================================

# DFG v0.18/19
export FunctorInferenceType, PackedInferenceType

@deprecate _evalType(pt::String) DFG.getTypeFromSerializationModule(pt)


export InMemDFGType
const InMemDFGType = DFG.LocalDFG{SolverParams}

##==============================================================================
## Deprecate code below before v0.28
##==============================================================================

function Base.convert(::Type{String}, 
                      obj::FluxModelsDistribution)
  #
  @error "Obsolete, FluxModelsSerialization should not return String for general cases of PackedSamplableBelief"
  # convert to packed type first
  packed = convert(PackedFluxModelsDistribution, obj)
  # FIXME, should not return String for general cases of PackedSamplableBelief 
  return JSON2.write(packed)
end


# import IncrementalInference: decodefg, loadjld

function veeCategorical(val::Categorical)
  @warn "veeCategorical is obsolete and being deprecated."
  val.p
end
function veeCategorical(val::Union{Nothing, Vector{Float64}})
  @warn "veeCategorical is obsolete and being deprecated."  
  val
end

function packmultihypo(fnc::CommonConvWrapper{T}) where {T<:AbstractFactor}
  @warn "packmultihypo is deprecated in favor of Vector only operations"
  fnc.hypotheses !== nothing ? string(fnc.hypotheses) : ""
end
function parsemultihypostr(str::AS) where {AS <: AbstractString}
  @warn "parsemultihypostr is deprecated in favor of Vector only operations"
  mhcat=nothing
  if length(str) > 0
    mhcat = convert(SamplableBelief, str)
  end
  return mhcat
end

# # FIXME DEPRECATE TO BETTER JSON with ._type field STANDARD
# function convert(::Type{<:PackedSamplableBelief}, obj::SamplableBelief)
#   # FIXME, prep for switch
#   packDistribution(obj)
  
#   # FIXME must use string, because unpacking templated e.g. PackedType{T} has problems, see DFG #668
#   string(obj)
# end


# New features towards standardizing distribution serialization
# # Assumes DFG/IIF serialized distributions have a `PackedType._type::String = "MyModule.MyPackedDistributionDensityType"`
# # also see DFG #590
# function convert( ::Type{String}, 
#                   obj::PackedSamplableBelief )
#   #
#   _typ = DFG.getTypeFromSerializationModule(obj._type)
# end

# convert(::Union{Type{<:SamplableBelief},Type{<:HeatmapGridDensity}},
#         obj::PackedHeatmapGridDensity) = unpackDistribution(obj)

# convert(::Union{Type{<:PackedSamplableBelief},Type{<:PackedHeatmapGridDensity}}, 
#         obj::HeatmapGridDensity ) = packDistribution(obj)
# #

# convert(::Union{Type{<:SamplableBelief},Type{<:LevelSetGridNormal}}, 
#         obj::PackedLevelSetGridNormal) = unpackDistribution(obj)

# convert(::Union{Type{<:PackedSamplableBelief},Type{<:PackedLevelSetGridNormal}}, 
#         obj::LevelSetGridNormal) = packDistribution(obj)


# TODO stop-gap string storage of Distrubtion types, should be upgraded to more efficient storage
function normalfromstring(str::AbstractString)
  meanstr = match(r"μ=[+-]?([0-9]*[.])?[0-9]+", str).match
  mean = split(meanstr, '=')[2]
  sigmastr = match(r"σ=[+-]?([0-9]*[.])?[0-9]+", str).match
  sigma = split(sigmastr, '=')[2]
  Normal{Float64}(parse(Float64,mean), parse(Float64,sigma))
end

function mvnormalfromstring(str::AbstractString)
  means = split(split(split(str, 'μ')[2],']')[1],'[')[end]
  mean = Float64[]
  for ms in split(means, ',')
    push!(mean, parse(Float64, ms))
  end
  sigs = split(split(split(str, 'Σ')[2],']')[1],'[')[end]
  sig = Float64[]
  for ms in split(sigs, ';')
    for m in split(ms, ' ')
      length(m) > 0 ? push!(sig, parse(Float64, m)) : nothing
    end
  end
  len = length(mean)
  sigm = reshape(sig, len,len)
  MvNormal(mean, sigm)
end

function categoricalfromstring(str::AbstractString)
  # pstr = match(r"p=\[", str).match
  psubs = split(str, '=')[end]
  psubs = split(psubs, '[')[end]
  psubsub = split(psubs, ']')[1]
  pw = split(psubsub, ',')
  p = parse.(Float64, pw)
  return Categorical(p ./ sum(p))
end



# NOTE SEE EXAMPLE IN src/Flux/FluxModelsSerialization.jl
function _extractDistributionJson(jsonstr::AbstractString, checkJson::AbstractVector{<:AbstractString})
  # Assume first word after split is the type
  mjs = findfirst(r"([a-zA-Z0-9._]+)\w", checkJson[2])
  maybemodule = split(checkJson[2][mjs], '.')
  # get dedicated Module or revert to Main
  packmodule = 1 < length(maybemodule) ? getfield(Main, Symbol(maybemodule[1])) : Main
  packtype = getfield(packmodule, Symbol(maybemodule[end]))
  packed = JSON2.read(jsonstr, packtype)
  # call the dedicated converter for this packed type using dispatch
  convert(SamplableBelief, packed)
end


function _legacyUnpackDistribution(str::Union{<:PackedSamplableBelief,<:AbstractString})
  # TODO improve use of multidispatch and packing of Distribution types
  # extractdistribution(str::AS) where {AS <: AbstractString}
  # TODO improve use of multidispatch and packing of Distribution types
  # TODO use startswith
  checkJson = split(str, r"PackedSamplableTypeJSON")
  if str == ""
    return nothing
  elseif length(checkJson) == 2
    # TODO this is the new direction for serializing (pack/unpack) of <:Samplable objects
    # NOTE uses intermediate consolidation keyword search pattern `SamplableTypeJSON`
    return _extractDistributionJson(str, checkJson)
  elseif occursin(r"_type", str) && occursin(r"ManifoldKernelDensity", str)
    return convert(ManifoldKernelDensity, str)
  elseif startswith(str, "DiagNormal")
    # Diags are internally squared, so only option here is to sqrt on input.
    return mvnormalfromstring(str)
  elseif startswith(str, "ZeroMeanDiagNormal")
    error("ZeroMeanDiagNormal not yet supported, deferring to full JSON serialization of all Distribution objects.")
  elseif occursin(r"FullNormal", str)
    return mvnormalfromstring(str)
  elseif (occursin(r"Normal", str) )# && !occursin(r"FullNormal", str))
    return normalfromstring(str)
  elseif occursin(r"Categorical", str)
    return categoricalfromstring(str)
  elseif occursin(r"DiscreteNonParametric", str)
    return categoricalfromstring(str)
  elseif occursin(r"KDE:", str)
    return convert(BallTreeDensity, str)
  elseif occursin(r"AliasingScalarSampler", str)
    return convert(AliasingScalarSampler, str)
  else
    error("Don't know how to extract distribution from str=$(str)")
  end
end

@deprecate convert(::Type{<:SamplableBelief}, str::Union{<:PackedSamplableBelief,<:AbstractString}) _legacyUnpackDistribution(str)

@deprecate HeatmapDensityRegular(w...;kw...) LevelSetGridNormal(w...;kw...)

@deprecate generateCanonicalFG_Kaess(w...;kw...) generateGraph_Kaess(w...;kw...)
@deprecate generateCanonicalFG_EuclidDistance(w...;kw...) generateGraph_EuclidDistance(w...;kw...)
@deprecate generateCanonicalFG_lineStep(w...;kw...) generateGraph_LineStep(w...;kw...)
@deprecate generateCanonicalFG_CaesarRing1D(w...;kw...) generateGraph_CaesarRing1D(w...;kw...)
@deprecate generateCanonicalFG_TestSymbolic(w...;kw...) generateGraph_TestSymbolic(w...;kw...)


##==============================================================================
## Deprecate code below before v0.27
##==============================================================================

convert(::Type{<:SamplableBelief}, obj::PackedUniform) = unpackDistribution(obj)


Base.convert(::Type{<:Tuple}, ::InstanceType{Manifolds.Euclidean{Tuple{N}, ℝ}} ) where N = tuple([:Euclid for i in 1:N]...)
Base.convert(::Type{<:Tuple}, ::InstanceType{Manifolds.Circle{ℝ}})  = error("#FIXME")#(:Circular,)
Base.convert(::Type{<:Tuple}, ::InstanceType{Manifolds.RealCircleGroup})  = (:Circular,)

@deprecate manikde!(pts::AbstractVector, bws::Vector{<:Real}, variableType::Union{InstanceType{<:InferenceVariable}, InstanceType{<:AbstractFactor}} ) manikde!(variableType,pts,bws)
@deprecate manikde!(pts::AbstractVector, varType::Union{InstanceType{<:InferenceVariable}, InstanceType{<:AbstractFactor}}) manikde!(varType, pts)

# function getSample(cf::CalcFactor{<:AbstractRelativeRoots})
#   M = getManifold(cf.factor)
#   if hasfield(typeof(cf.factor), :Z)
#     return sampleTangent(M, cf.factor.Z)
#   else
#     error("""Factor $(typeof(cf.factor)) does not have a field `Z`, to use the default `getSample` method, use `Z` for the measurement. 
#               Alternatively, provide a `getSample` method. See IIF issue #1441 and Custom Factors in the Caesar documentation.""")
#   end
# end

# function getSample(cf::CalcFactor{<:CircularCircular})
#   # return (sampleTangent(getManifold(cf.factor), cf.factor.Z), )
#   # FIXME workaround for issue #TBD with manifolds CircularGroup
#   return sampleTangent(getManifold(cf.factor), cf.factor.Z, [0.0])
# end

# function getSample(cf::CalcFactor{<:LinearRelative})
#   # _samplemakevec(z::Real) = [z;]
#   # _samplemakevec(z::AbstractVector{<:Real}) = z
#   return sampleTangent(getManifold(cf.factor), cf.factor.Z)
# end

# function getSample(cf::CalcFactor{<:EuclidDistance})
#   rand(cf.factor.Z, 1)
# end

@deprecate CalcFactor(x1,x2,x3,x4,x5,x6) CalcFactor(x1,x2,x3,x4,x5,x6, true)

@deprecate ensureAllInitialized!(w...;kw...) initAll!(w...;kw...)

# """
#     $SIGNATURES

# Sample points from a regular grid scalar field level set.

# Notes
# - modifies first argument roi to save on memory allocations
# - user should duplicate if with a deepcopy if needed
# """
# function sampleLevelSetGaussian!( roi::AbstractMatrix{<:Real},
#                                   sigma::Real,
#                                   x_grid::AbstractVector{<:Real}, 
#                                   y_grid::AbstractVector{<:Real};
#                                   sigma_scale::Real=3  )
#   #
#   # make Gaussian
#   roi .^= 2
#   roi .*= 0.5/(sigma^2)
#     # truncate at sigma_scale*sigma
#     _roi = thres .- roi

#   sampleHeatmap(roi, x_grid, y_grid, sigma_scale^2)
# end


# @deprecate getLevelSetSigma(data::AbstractMatrix{<:Real}, level::Real, w...; kw...) sampleLevelSetGaussian!(data.-level, w...; kw...)

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



#
