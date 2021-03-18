
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
# Sphere1Sphere1()


@deprecate Sphere1Sphere1(w...;kw...) CircularCircular(w...;kw...)

getSample(cf::CalcFactor{<:CircularCircular}, N::Int=1) = (reshape(rand(cf.factor.Z,N),:,N), )


function (cf::CalcFactor{<:CircularCircular})(meas,
                                            wxi,
                                            wxj  ) # where {M<:FactorMetadata,P<:Tuple,X<:AbstractVector}
  #
  wXjhat = addtheta(wxi[1], meas[1])
  res = difftheta(wxj[1], wXjhat)  # jXjhat =
  return res
end


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

@deprecate PriorSphere1(w...;kw...) PriorCircular(w...;kw...)

PriorCircular(::UniformScaling) = PriorCircular(Normal())


function getSample(cf::CalcFactor{<:PriorCircular}, N::Int=1)
  return (reshape(rand(cf.factor.Z,N),:,N), )
end





"""
$(TYPEDEF)

Serialized object for storing PriorCircular.
"""
mutable struct PackedPriorCircular  <: IncrementalInference.PackedInferenceType
  datastr::String
  PackedPriorCircular() = new()
  PackedPriorCircular(x::String) = new(x)
end
function convert(::Type{PackedPriorCircular}, d::PriorCircular)
  return PackedPriorCircular(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{PriorCircular}, d::PackedPriorCircular)
  distr = convert(SamplableBelief, d.datastr)
  return PriorCircular{typeof(distr)}(distr)
end

@deprecate PackedPriorSphere1(w...;kw...) PackedPriorCircular(w...;kw...)

# --------------------------------------------





"""
$(TYPEDEF)

Serialized object for storing CircularCircular.
"""
mutable struct PackedCircularCircular  <: IncrementalInference.PackedInferenceType
  datastr::String
  PackedCircularCircular() = new()
  PackedCircularCircular(x::String) = new(x)
end
function convert(::Type{CircularCircular}, d::PackedCircularCircular)
  return CircularCircular(convert(SamplableBelief, d.datastr))
end
function convert(::Type{PackedCircularCircular}, d::CircularCircular)
  return PackedCircularCircular(convert(PackedSamplableBelief, d.Z))
end


@deprecate PackedSphere1Sphere1(w...;kw...) PackedCircularCircular(w...;kw...)


# --------------------------------------------
