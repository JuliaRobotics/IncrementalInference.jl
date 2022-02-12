
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

@deprecate kde!(em::TreeBelief) manikde!(em)

# DFG v0.18/19
export FunctorInferenceType, PackedInferenceType

@deprecate _evalType(pt::String) DFG.getTypeFromSerializationModule(pt)

# LightDFG will be replaced by GraphsDFG
export LightDFG
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




#
