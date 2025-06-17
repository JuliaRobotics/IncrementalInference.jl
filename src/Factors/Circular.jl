
export CircularCircular, PriorCircular, PackedCircularCircular, PackedPriorCircular

"""
$(TYPEDEF)

Factor between two Sphere1 variables.

Related

[`Sphere1`](@ref), [`PriorSphere1`](@ref), [`Polar`](@ref), [`ContinuousEuclid`](@ref)
"""
mutable struct CircularCircular{T <: SamplableBelief} <: AbstractManifoldMinimize
  Z::T
  # Sphere1Sphere1(z::T=Normal()) where {T <: SamplableBelief} = new{T}(z)
end

const Sphere1Sphere1 = CircularCircular

CircularCircular(::UniformScaling) = CircularCircular(Normal())

DFG.getManifold(::CircularCircular) = RealCircleGroup()

function (cf::CalcFactor{<:CircularCircular})(X, p, q)
  #
  M = getManifold(cf)
  return distanceTangent2Point(M, X, p, q)
end

function getSample(cf::CalcFactor{<:CircularCircular})
  # FIXME workaround for issue with manifolds CircularGroup, 
  return [rand(cf.factor.Z)]
end

function Base.convert(::Type{<:MB.AbstractManifold}, ::InstanceType{CircularCircular})
  return Manifolds.RealCircleGroup()
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
DFG.@defFactorType PriorCircular AbstractPrior Manifolds.RealCircleGroup()

# mutable struct PriorCircular{T <: SamplableBelief} <: AbstractPrior
#   Z::T
# end

PriorCircular(::UniformScaling) = PriorCircular(Normal())

# DFG.getManifold(::PriorCircular) = RealCircleGroup()

function getSample(cf::CalcFactor{<:PriorCircular})
  # FIXME workaround for issue #TBD with manifolds CircularGroup, 
  # JuliaManifolds/Manifolds.jl#415
  # no method similar(::Float64, ::Type{Float64})
  return samplePoint(cf.manifold, cf.factor.Z, [0.0])
  # return [Manifolds.sym_rem(rand(cf.factor.Z))]
end

function (cf::CalcFactor{<:PriorCircular})(m, p)
  M = getManifold(cf)
  Xc = vee(M, p, log(M, p, m))
  return Xc
end

function Base.convert(::Type{<:MB.AbstractManifold}, ::InstanceType{PriorCircular})
  return Manifolds.RealCircleGroup()
end

# """
# $(TYPEDEF)

# Serialized object for storing PriorCircular.
# """
# Base.@kwdef struct PackedPriorCircular <: AbstractPackedFactor
#   Z::PackedSamplableBelief
# end

# function convert(::Type{PackedPriorCircular}, d::PriorCircular)
#   return PackedPriorCircular(convert(PackedSamplableBelief, d.Z))
# end
# function convert(::Type{PriorCircular}, d::PackedPriorCircular)
#   distr = convert(SamplableBelief, d.Z)
#   return PriorCircular{typeof(distr)}(distr)
# end

# --------------------------------------------

"""
$(TYPEDEF)

Serialized object for storing CircularCircular.
"""
Base.@kwdef struct PackedCircularCircular <: AbstractPackedFactor
  Z::PackedSamplableBelief
end
function convert(::Type{CircularCircular}, d::PackedCircularCircular)
  return CircularCircular(convert(SamplableBelief, d.Z))
end
function convert(::Type{PackedCircularCircular}, d::CircularCircular)
  return PackedCircularCircular(convert(PackedSamplableBelief, d.Z))
end

# --------------------------------------------
