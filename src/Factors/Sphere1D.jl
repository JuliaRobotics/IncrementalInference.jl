
export Sphere1Sphere1, PriorSphere1, PackedSphere1Sphere1, PackedPriorSphere1


"""
$(TYPEDEF)
"""
mutable struct Sphere1Sphere1{T} <: IncrementalInference.FunctorPairwise where {T <: SamplableBelief}
  Z::T
  Sphere1Sphere1{T}() where {T <: SamplableBelief} = new{T}()
  Sphere1Sphere1{T}(z1::T) where {T <: SamplableBelief} = new{T}(z1)
end
Sphere1Sphere1(z::T) where {T <: SamplableBelief} = Sphere1Sphere1{T}(z)

getSample(s::Sphere1Sphere1{<: SamplableBelief}, N::Int=1) = (rand(s.Z, N), )

function (s::Sphere1Sphere1{<: SamplableBelief})(res::Vector{Float64},
                                                 userdata::FactorMetadata,
                                                 idx::Int,
                                                 meas::Tuple,
                                                 wxi::Array{Float64,2},
                                                 wxj::Array{Float64,2}  )
  #
  wXjhat = addtheta(wxi[1,idx], meas[1][1,idx])
  res[1] = difftheta(wxj[1,idx], wXjhat)  # jXjhat =
  nothing
end


"""
$(TYPEDEF)

Introduce direct observations on all dimensions of a Sphere1 variable:

Example:
--------
```julia
PriorSphere1( MvNormal([10; 10; pi/6.0], Matrix(Diagonal([0.1;0.1;0.05].^2))) )
```
"""
mutable struct PriorSphere1{T} <: IncrementalInference.FunctorSingleton  where {T <: IncrementalInference.SamplableBelief}
    Z::T
    PriorSphere1{T}() where T = new{T}()
    PriorSphere1{T}(x::T) where {T <: IncrementalInference.SamplableBelief}  = new{T}(x)
end
PriorSphere1(x::T) where {T <: IncrementalInference.SamplableBelief} = PriorSphere1{T}(x)
function PriorSphere1(mu::Array{Float64}, cov::Array{Float64,2}, W::Vector{Float64})
  @warn "PriorSphere1(mu,cov,W) is deprecated in favor of PriorSphere1(T(...)) -- use for example PriorSphere1(MvNormal(mu, cov))"
  PriorSphere1(MvNormal(mu[:], cov))
end
function getSample(p2::PriorSphere1, N::Int=1)
  return (rand(p2.Z,N), )
end





"""
$(TYPEDEF)

Serialized object for storing PriorSphere1.
"""
mutable struct PackedPriorSphere1  <: IncrementalInference.PackedInferenceType
  datastr::String
  PackedPriorSphere1() = new()
  PackedPriorSphere1(x::String) = new(x)
end
function convert(::Type{PackedPriorSphere1}, d::PriorSphere1)
  return PackedPriorSphere1(string(d.Z))
end
function convert(::Type{PriorSphere1}, d::PackedPriorSphere1)
  distr = extractdistribution(d.datastr)
  return PriorSphere1{typeof(distr)}(distr)
end



# --------------------------------------------





"""
$(TYPEDEF)

Serialized object for storing Sphere1Sphere1.
"""
mutable struct PackedSphere1Sphere1  <: IncrementalInference.PackedInferenceType
  datastr::String
  PackedSphere1Sphere1() = new()
  PackedSphere1Sphere1(x::String) = new(x)
end
function convert(::Type{Sphere1Sphere1}, d::PackedSphere1Sphere1)
  return Sphere1Sphere1(extractdistribution(d.datastr))
end
function convert(::Type{PackedSphere1Sphere1}, d::Sphere1Sphere1)
  return PackedSphere1Sphere1(string(d.z))
end





# --------------------------------------------
