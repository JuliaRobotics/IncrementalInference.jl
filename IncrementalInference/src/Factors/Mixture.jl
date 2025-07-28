
_defaultNamesMixtures(N::Int) = ((Symbol[Symbol("c$i") for i = 1:N])...,)

"""
$(TYPEDEF)

A `Mixture` object for use with either a `<: AbstractPrior` or `<: AbstractRelative`.

Notes
- The internal data representation is a `::NamedTuple`, which allows total type-stability for all component types.
- Various construction helpers can accept a variety of inputs, including `<: AbstractArray` and `Tuple`.
- `N` is the number of components used to make the mixture, so two bumps from two Normal components means `N=2`.

DevNotes
- FIXME swap API order so Mixture of distibutions works like a distribtion, see Caesar.jl #808
  - Should not have field mechanics.
- TODO on sampling see #1099 and #1094 and #1069 


Example
```juila
# prior factor
msp = Mixture(Prior, 
              [Normal(0,0.1), Uniform(-pi/1,pi/2)],
              [0.5;0.5])

addFactor!(fg, [:head], msp, tags=[:MAGNETOMETER;])

# Or relative
mlr = Mixture(LinearRelative, 
              (correlator=AliasingScalarSampler(...), naive=Normal(0.5,5), lucky=Uniform(0,10)),
              [0.5;0.4;0.1])

addFactor!(fg, [:x0;:x1], mlr)
```
"""
struct Mixture{N, F <: AbstractFactor, S, T <: Tuple} <: AbstractFactor
  """ factor mechanics """
  mechanics::F
  components::NamedTuple{S, T}
  diversity::Distributions.Categorical
  """ dimension of factor, so range measurement would be dims=1 """
  dims::Int
  labels::Vector{Int}
end

function Mixture(
  f::Type{F},
  z::NamedTuple{S, T},
  c::Distributions.DiscreteNonParametric,
) where {F <: AbstractFactor, S, T}
  return Mixture{length(z), F, S, T}(
    f(LinearAlgebra.I),
    z,
    c,
    size(rand(z[1], 1), 1),
    zeros(Int, 0),
  )
end
function Mixture(
  f::F,
  z::NamedTuple{S, T},
  c::Distributions.DiscreteNonParametric,
) where {F <: AbstractFactor, S, T}
  return Mixture{length(z), F, S, T}(f, z, c, size(rand(z[1], 1), 1), zeros(Int, 0))
end
function Mixture(
  f::Union{F, Type{F}},
  z::NamedTuple{S, T},
  c::AbstractVector{<:Real},
) where {F <: AbstractFactor, S, T}
  return Mixture(f, z, Categorical([c...]))
end
function Mixture(
  f::Union{F, Type{F}},
  z::NamedTuple{S, T},
  c::NTuple{N, <:Real},
) where {N, F <: AbstractFactor, S, T}
  return Mixture(f, z, [c...])
end
function Mixture(
  f::Union{F, Type{F}},
  z::Tuple,
  c::Union{
    <:Distributions.DiscreteNonParametric,
    <:AbstractVector{<:Real},
    <:NTuple{N, <:Real},
  },
) where {F <: AbstractFactor, N}
  return Mixture(f, NamedTuple{_defaultNamesMixtures(length(z))}(z), c)
end
function Mixture(
  f::Union{F, Type{F}},
  z::AbstractVector{<:SamplableBelief},
  c::Union{
    <:Distributions.DiscreteNonParametric,
    <:AbstractVector{<:Real},
    <:NTuple{N, <:Real},
  },
) where {F <: AbstractFactor, N}
  return Mixture(f, (z...,), c)
end

function Base.resize!(mp::Mixture, s::Int)
  return resize!(mp.labels, s)
end

_lengthOrNothing(val) = length(val)
_lengthOrNothing(val::Nothing) = 0

getManifold(m::Mixture) = getManifold(m.mechanics)

# TODO make in-place memory version
function sampleFactor(cf::CalcFactor{<:Mixture}, N::Int = 1)
  #
  # TODO consolidate #927, case if mechanics has a special sampler
  # TODO slight bit of waste in computation, but easiest way to ensure special tricks in s.mechanics::F are included
  ## example case is old FluxModelsPose2Pose2 requiring velocity
  # FIXME better consolidation of when to pass down .mechanics, also see #1099 and #1094 and #1069

  cf_ = CalcFactorNormSq(
    cf.factor.mechanics,
    0,
    cf._legacyParams,
    cf._allowThreads,
    cf.cache,
    cf.fullvariables,
    cf.solvefor,
    cf.manifold,
    cf.measurement,
    nothing,
  )
  smpls = [getSample(cf_) for _ = 1:N]
  # smpls = Array{Float64,2}(undef,s.dims,N)
  #out memory should be right size first
  length(cf.factor.labels) != N ? resize!(cf.factor.labels, N) : nothing
  cf.factor.labels .= rand(cf.factor.diversity, N)
  M = cf.manifold
  
  # mixture needs to be refactored so let's make it worse :-)
  if cf.factor.mechanics isa AbstractPrior
    samplef = samplePoint
  elseif cf.factor.mechanics isa AbstractRelative
    samplef = sampleTangent
  end

  for i = 1:N
    mixComponent = cf.factor.components[cf.factor.labels[i]]
    # measurements relate to the factor's manifold (either tangent vector or manifold point)
    setPointsMani!(smpls,  samplef(M, mixComponent), i)
  end

  # TODO only does first element of meas::Tuple at this stage, see #1099
  return smpls
end

function DistributedFactorGraphs.isPrior(::Type{Mixture{N, F, S, T}}) where {N, F, S, T}
  return F <: AbstractPriorObservation
end

"""
$(TYPEDEF)

Serialization type for `Mixture`.
"""
Base.@kwdef mutable struct PackedMixture <: AbstractPackedFactor
  N::Int
  # store the packed type for later unpacking
  F_::String
  S::Vector{String}
  components::Vector{PackedBelief}
  diversity::PackedBelief
end

function convert(::Type{<:PackedMixture}, obj::Mixture{N, F, S, T}) where {N, F, S, T}
  allcomp = PackedBelief[]
  for val in obj.components
    dtr_ = convert(PackedBelief, val)
    # FIXME ON FIRE, likely to be difficult for non-standard "Samplable" types -- e.g. Flux models in RoME
    push!(allcomp, dtr_)
  end
  if hasmethod(pack, (typeof(obj.mechanics),))
    pm = pack(obj.mechanics)
  else
    @warn("No pack method for mechanics type $(typeof(obj.mechanics)), using deprecated convert instead.")
    pm = convert(DFG.convertPackedType(obj.mechanics), obj.mechanics)
  end
  sT = string(typeof(pm))
  dvst = convert(PackedBelief, obj.diversity)
  return PackedMixture(N, sT, string.(collect(S)), allcomp, dvst)
end

function convert(::Type{<:Mixture}, obj::PackedMixture)
  N = obj.N
  F1 = getfield(Main, Symbol(obj.F_))
  S = (Symbol.(obj.S)...,)
  F2 = DFG.convertStructType(F1)

  components = map(c -> convert(SamplableBelief, c), obj.components)
  diversity = convert(SamplableBelief, obj.diversity)
  # tupcomp = (components...,)
  ntup = NamedTuple{S}(components) # ,typeof(tupcomp)
  return Mixture(F2, ntup, diversity)
end

#
