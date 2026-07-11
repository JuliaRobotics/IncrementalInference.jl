"""
$(TYPEDEF)

A `Mixture` object for use with either a `<: AbstractPriorObservation` or `<: AbstractRelativeObservation`.

Notes
- The internal data representation is a `::NamedTuple`, which allows total type-stability for all component types.
- Various construction helpers can accept a variety of inputs, including `<: AbstractArray` and `Tuple`.
- `N` is the number of components used to make the mixture, so two bumps from two Normal components means `N=2`.

DevNotes
- TODO on sampling see #1099 and #1094 and #1069 

Example
```julia
# prior factor
msp = Prior(Mixture([Normal(0,0.1), Uniform(-pi/1,pi/2)], [0.5;0.5])

addFactor!(fg, [:head], msp, tags=[:MAGNETOMETER;])

# Or relative
mlr = LinearRelative(Mixture(
              (correlator=AliasingScalarSampler(...), naive=Normal(0.5,5), lucky=Uniform(0,10)),
              [0.5;0.4;0.1]
      )

addFactor!(fg, [:x0;:x1], mlr)
```
"""
# we can do <: Sampleable{VF, VS}, but then all components must implement Sampleable
struct Mixture{C <: NamedTuple, P <: Categorical}
    components::C
    prior::P
end

Base.length(d::Mixture) = length(d.components[1])
#FIXME next not working
# Base.rand(rng::AbstractRNG, d::Mixture) = rand(rng, d.components[rand(rng, d.prior)])
Base.rand(d::Mixture) = rand(d.components[rand(d.prior)])
getDimension(p::Mixture) = length(p)

Mixture(z::NamedTuple) = Mixture(z, Categorical(length(z)))
Mixture(z::NamedTuple, c::Union{<:Tuple, <:AbstractVector}) = Mixture(z, Categorical(c))

function Mixture(
    z::Union{<:Tuple, <:AbstractVector},
    args...,
)
    cnames = tuple((Symbol("c", i) for i = 1:length(z))...)
    return Mixture(NamedTuple{cnames}(z), args...)
end

function Mixture(f::Type{<:AbstractObservation}, args...)
    return f(Mixture(args...))
end

DFG.@tags struct PackedMixture
  components::NamedTuple & (choosetype = obj->NamedTuple{Tuple(Symbol.(keys(obj))), Tuple{DFG.resolveType.(values(obj))...}}, )
  prior & DFG.@packed
end

function DFG.pack(m::Mixture)
  PackedMixture(map(DFG.Packed, m.components), m.prior)
end

function DFG.unpack(pm::PackedMixture)
  Mixture(
    map(unpack, pm.components),
    pm.prior
  )
end
