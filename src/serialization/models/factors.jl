export
    Packed_Factor

mutable struct Packed_Factor
    measurement::Vector{Dict{String, Any}}
    additionalData::Vector{Dict{String, Any}}
    factorType::String
end
