
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
## Deprecate code below before v0.31
##==============================================================================


@deprecate initManual!(w...;kw...) initVariable!(w...;kw...)



##==============================================================================
## Deprecate code below before v0.29
##==============================================================================

# function manikde!(vartype::Union{InstanceType{<:InferenceVariable}, InstanceType{<:AbstractFactor}},
#                   pts::AbstractVector{P} ) where {P <: Union{<:AbstractArray,<:Number,<:ProductRepr,<:Manifolds.ArrayPartition} }
#   #
#   M = getManifold(vartype)
#   infoPerCoord=ones(manifold_dimension(M))
#   return AMP.manikde!(M, pts; infoPerCoord)
# end

# @deprecate manikde!(varType::Union{InstanceType{<:InferenceVariable}, InstanceType{<:AbstractFactor}}, pts::AbstractVector{P}, bw::Vector{<:Real} ) where {P <: Union{<:AbstractArray,<:Number,<:ProductRepr,<:Manifolds.ArrayPartition} } manikde!(varType, pts; bw)

# @deprecate kde!(em::TreeBelief) manikde!(em)

# # DFG v0.18/19
# export FunctorInferenceType, PackedInferenceType

# @deprecate _evalType(pt::String) DFG.getTypeFromSerializationModule(pt)

# # LightDFG will be replaced by GraphsDFG
# export LightDFG
# export InMemDFGType
# const InMemDFGType = DFG.LocalDFG{SolverParams}
