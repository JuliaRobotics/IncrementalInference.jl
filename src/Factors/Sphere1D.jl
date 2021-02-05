
export Sphere1Sphere1, PriorSphere1, PackedSphere1Sphere1, PackedPriorSphere1


"""
$(TYPEDEF)

Factor between two Sphere1 variables.

Related

[`Sphere1`](@ref), [`PriorSphere1`](@ref), [`Polar`](@ref), [`ContinuousEuclid`](@ref)
"""
mutable struct Sphere1Sphere1{T<: SamplableBelief} <: AbstractRelativeRoots
  Z::T
  Sphere1Sphere1(z::T=Normal()) where {T <: SamplableBelief} = new{T}(z)
end
Sphere1Sphere1(::UniformScaling) = Sphere1Sphere1()


getSample(cf::CalcFactor{<:Sphere1Sphere1}, N::Int=1) = (reshape(rand(cf.factor.Z,N),:,N), )


function (cf::CalcFactor{<:Sphere1Sphere1})(meas,
                                            wxi,
                                            wxj  ) # where {M<:FactorMetadata,P<:Tuple,X<:AbstractVector}
  #
  wXjhat = addtheta(wxi[1], meas[1])
  res = difftheta(wxj[1], wXjhat)  # jXjhat =
  return res
end


"""
$(TYPEDEF)

Introduce direct observations on all dimensions of a Sphere1 variable:

Example:
--------
```julia
PriorSphere1( MvNormal([10; 10; pi/6.0], diagm([0.1;0.1;0.05].^2)) )
```

Related

[`Sphere1`](@ref), [`Prior`](@ref), [`PartialPrior`](@ref)
"""
mutable struct PriorSphere1{T<: SamplableBelief} <: AbstractPrior
  Z::T
end


# PriorSphere1(x::T) where {T <: IncrementalInference.SamplableBelief} = PriorSphere1{T}(x)
PriorSphere1(::UniformScaling) = PriorSphere1(Normal())


function getSample(cf::CalcFactor{<:PriorSphere1}, N::Int=1)
  return (reshape(rand(cf.factor.Z,N),:,N), )
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
  return PackedPriorSphere1(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{PriorSphere1}, d::PackedPriorSphere1)
  distr = convert(SamplableBelief, d.datastr)
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
  return Sphere1Sphere1(convert(SamplableBelief, d.datastr))
end
function convert(::Type{PackedSphere1Sphere1}, d::Sphere1Sphere1)
  return PackedSphere1Sphere1(convert(PackedSamplableBelief, d.Z))
end





# --------------------------------------------
