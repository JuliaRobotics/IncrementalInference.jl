# basic agnostic factors.  See RoME.jl for more comprehensive factors

"""
$(TYPEDEF)
"""
mutable struct Ranged <: AbstractRelativeRoots
    Zij::Array{Float64,1}
    Cov::Array{Float64,1}
    W::Array{Float64,1}
    Ranged() = new()
    Ranged(x...) = new(x[1], x[2], x[3])
end
"""
$(TYPEDEF)
"""
mutable struct PackedRanged <: PackedInferenceType
    Zij::Array{Float64,1}
    Cov::Array{Float64,1}
    W::Array{Float64,1}
    PackedRanged() = new()
    PackedRanged(x...) = new(x[1], x[2], x[3])
end
function convert(::Type{Ranged}, r::PackedRanged)
  return Ranged(r.Zij, r.Cov, r.W)
end
function convert(::Type{PackedRanged}, r::Ranged)
  return PackedRanged(r.Zij, r.Cov, r.W)
end
function (ra::Ranged)(res::AbstractVector{<:Real},
    userdata::FactorMetadata,
    idx::Int,
    meas::Tuple{<:AbstractArray{<:Real,2}},
    p1::AbstractArray{<:Real},
    l1::AbstractArray{<:Real})

  res[1] = meas[1][1,idx] - abs(l1[1,idx] - p1[1,idx])
  nothing
end
function getSample(ra::Ranged, N::Int=1)
  ret = zeros(1,N)
  for i in 1:N
    ret[1,i] = ra.Cov[1]*randn()+ra.Zij[1]
  end
  # rand(Distributions.Normal(odo.Zij[1],odo.Cov[1]), N)'
  return (ret,)
end


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
