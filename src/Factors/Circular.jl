
export CircularCircular, PriorCircular, PackedCircularCircular, PackedPriorCircular


"""
$(TYPEDEF)

Factor between two Sphere1 variables.

Related

[`Sphere1`](@ref), [`PriorSphere1`](@ref), [`Polar`](@ref), [`ContinuousEuclid`](@ref)
"""
mutable struct CircularCircular{T<: SamplableBelief} <: AbstractRelativeRoots
  Z::T
  # Sphere1Sphere1(z::T=Normal()) where {T <: SamplableBelief} = new{T}(z)
end

const Sphere1Sphere1 = CircularCircular

CircularCircular(::UniformScaling) = CircularCircular(Normal())


function getSample(cf::CalcFactor{<:CircularCircular}, N::Int=1)
  (randToPoints(cf.factor.Z, N), )
end

function (cf::CalcFactor{<:CircularCircular})(meas,
                                            wxi,
                                            wxj  ) # where {M<:FactorMetadata,P<:Tuple,X<:AbstractVector}
  #
  wXjhat = addtheta(wxi[1], meas[1])
  res = difftheta(wxj[1], wXjhat)  # jXjhat =
  return res
end


Base.convert(::Type{<:MB.AbstractManifold}, ::InstanceType{CircularCircular}) = Manifolds.Circle()


"""
$(TYPEDEF)

Introduce direct observations on all dimensions of a Circular variable:

Example:
--------
```julia
PriorCircular( MvNormal([10; 10; pi/6.0], diagm([0.1;0.1;0.05].^2)) )
```

Related

[`Circular`](@ref), [`Prior`](@ref), [`PartialPrior`](@ref)
"""
mutable struct PriorCircular{T<: SamplableBelief} <: AbstractPrior
  Z::T
end

PriorCircular(::UniformScaling) = PriorCircular(Normal())


function getSample(cf::CalcFactor{<:PriorCircular}, N::Int=1)
  return (randToPoints(cf.factor.Z,N),)
end


Base.convert(::Type{<:MB.AbstractManifold}, ::InstanceType{PriorCircular}) = Manifolds.Circle()



"""
$(TYPEDEF)

Serialized object for storing PriorCircular.
"""
mutable struct PackedPriorCircular  <: IncrementalInference.PackedInferenceType
  datastr::String
  # PackedPriorCircular() = new()
  # PackedPriorCircular(x::String) = new(x)
end
function convert(::Type{PackedPriorCircular}, d::PriorCircular)
  return PackedPriorCircular(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{PriorCircular}, d::PackedPriorCircular)
  distr = convert(SamplableBelief, d.datastr)
  return PriorCircular{typeof(distr)}(distr)
end


# --------------------------------------------





"""
$(TYPEDEF)

Serialized object for storing CircularCircular.
"""
mutable struct PackedCircularCircular  <: IncrementalInference.PackedInferenceType
  datastr::String
  # PackedCircularCircular() = new()
  # PackedCircularCircular(x::String) = new(x)
end
function convert(::Type{CircularCircular}, d::PackedCircularCircular)
  return CircularCircular(convert(SamplableBelief, d.datastr))
end
function convert(::Type{PackedCircularCircular}, d::CircularCircular)
  return PackedCircularCircular(convert(PackedSamplableBelief, d.Z))
end



# --------------------------------------------
