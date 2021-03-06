# agnostic factors for internal use


"""
$(TYPEDEF)
"""
mutable struct GenericMarginal <: AbstractRelativeRoots
    Zij::Array{Float64,1}
    Cov::Array{Float64,1}
    W::Array{Float64,1}
    GenericMarginal() = new()
    GenericMarginal(a,b,c) = new(a,b,c)
end
"""
$(TYPEDEF)
"""
mutable struct PackedGenericMarginal <: PackedInferenceType
    Zij::Array{Float64,1}
    Cov::Array{Float64,1}
    W::Array{Float64,1}
    PackedGenericMarginal() = new()
    PackedGenericMarginal(a,b,c) = new(a,b,c)
end
function convert(::Type{PackedGenericMarginal}, d::GenericMarginal)
  return PackedGenericMarginal(d.Zij, d.Cov, d.W)
end
function convert(::Type{GenericMarginal}, d::PackedGenericMarginal)
  return GenericMarginal(d.Zij, d.Cov, d.W)
end

# ------------------------------------------------------------
