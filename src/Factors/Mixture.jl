
export Mixture, PackedMixture


_defaultNamesMixtures(N::Int) = ((Symbol[Symbol("c$i") for i in 1:N])...,)


"""
$(TYPEDEF)

A `Mixture` object for use with either a `<: AbstractPrior` or `<: AbstractRelative`.

Notes
- The internal data representation is a `::NamedTuple`, which allows total type-stability for all component types.
- Various construction helpers can accept a variety of inputs, including `<: AbstractArray` and `Tuple`.

DevNotes
- TODO on sampling see #1099 and #1094 and #1069 

Example
```juila
# prior factor
msp = Mixture(PriorSphere1, 
              [model=Normal(0,0.1), Uniform(-pi/1,pi/2)],
              [0.5;0.5])

addFactor!(fg, [:head], msp, tags=[:MAGNETOMETER;])

# Or relative
mlr = Mixture(LinearRelative, 
              (correlator=AliasingScalarSampler(...), naive=Normal(0.5,5), lucky=Uniform(0,10)),
              [0.5;0.4;0.1])

addFactor!(fg, [:x0;:x1], mlr)
```
"""
struct Mixture{N, F<:FunctorInferenceType, S, T<:Tuple} <: FunctorInferenceType
  mechanics::F
  components::NamedTuple{S,T}
  diversity::Distributions.Categorical
  dims::Int
  labels::Vector{Int}
end


Mixture(f::Type{F},
        z::NamedTuple{S,T}, 
        c::Distributions.DiscreteNonParametric ) where {F<:FunctorInferenceType, S, T} = Mixture{length(z),F,S,T}(f(LinearAlgebra.I), z, c, size( rand(z[1],1), 1), zeros(Int, 0))
Mixture(f::F,
        z::NamedTuple{S,T}, 
        c::Distributions.DiscreteNonParametric ) where {F<:FunctorInferenceType, S, T} = Mixture{length(z),F,S,T}(f, z, c, size( rand(z[1],1), 1), zeros(Int, 0))
Mixture(f::Union{F,Type{F}},z::NamedTuple{S,T}, 
        c::AbstractVector{<:Real}) where {F<:FunctorInferenceType,S,T} = Mixture(f, z, Categorical([c...]) )
Mixture(f::Union{F,Type{F}},
        z::NamedTuple{S,T}, 
        c::NTuple{N,<:Real}) where {N,F<:FunctorInferenceType,S,T} = Mixture(f, z, [c...] )
Mixture(f::Union{F,Type{F}},
        z::Tuple, 
        c::Union{<:Distributions.DiscreteNonParametric, <:AbstractVector{<:Real}, <:NTuple{N,<:Real}} ) where {F<:FunctorInferenceType, N} = Mixture(f,NamedTuple{_defaultNamesMixtures(length(z))}(z), c )
Mixture(f::Union{F,Type{F}},
        z::AbstractVector{<:SamplableBelief}, 
        c::Union{<:Distributions.DiscreteNonParametric, <:AbstractVector{<:Real}, <:NTuple{N,<:Real}} ) where {F <: FunctorInferenceType, N} = Mixture(f,(z...,), c )


function Base.resize!(mp::Mixture, s::Int)
  resize!(mp.labels, s)
end

# TODO make in-place memory version
function getSample( cf::CalcFactor{<:Mixture}, 
                    N::Int=1  )
  #
  # TODO consolidate #927, case if mechanics has a special sampler
  # TODO slight bit of waste in computation, but easiest way to ensure special tricks in s.mechanics::F are included
  ## example case is old FluxModelsPose2Pose2 requiring velocity
  # FIXME better consolidation of when to pass down .mechanics, also see #1099 and #1094 and #1069
  cf_ = CalcFactor( cf.factor.mechanics, cf.metadata, 0, length(cf._legacyMeas), cf._legacyMeas, cf._legacyParams)
  smpls = getSample(cf_, N)
    # smpls = Array{Float64,2}(undef,s.dims,N)
  #out memory should be right size first
  length(cf.factor.labels) != N ? resize!(cf.factor.labels, N) : nothing
  cf.factor.labels .= rand(cf.factor.diversity, N)
  for i in 1:N
    mixComponent = cf.factor.components[cf.factor.labels[i]]
    # measurements relate to the factor's manifold (either tangent vector or manifold point)
    _setPointsMani!(smpls[1][i], randToPoints(mixComponent,1)[1])
  end

  # TODO only does first element of meas::Tuple at this stage, see #1099
  smpls
end




function DistributedFactorGraphs.isPrior(::Mixture{N, F, S, T}) where {N,F,S,T}
  return F <: AbstractPrior
end

"""
$(TYPEDEF)

Serialization type for `Mixture`.
"""
mutable struct PackedMixture <: PackedInferenceType
  N::Int
  # store the packed type for later unpacking
  F_::String 
  S::Vector{String}
  components::Vector{String}
  diversity::String
end


function convert(::Type{PackedMixture}, obj::Mixture{N,F,S,T}) where {N,F,S,T}
  allcomp = String[]
  for val in obj.components
    # TODO, likely to be difficult for non-standard "Samplable" types -- e.g. Flux models in RoME
    push!(allcomp, convert(PackedSamplableBelief, val))
  end
  # pm = DFG.convertPackedType(obj.mechanics)
  pm = convert(DFG.convertPackedType(obj.mechanics), obj.mechanics)
  sT = string(typeof(pm))
  PackedMixture( N, sT, string.(collect(S)), allcomp, convert(PackedSamplableBelief, obj.diversity) )
end
function convert(::Type{Mixture}, obj::PackedMixture)
  N = obj.N
  F1 = getfield(Main, Symbol(obj.F_))
  S = (Symbol.(obj.S)...,)
  F2 = DFG.convertStructType(F1)
  components = convert.(SamplableBelief, obj.components)
  diversity = convert(SamplableBelief, obj.diversity)
  tupcomp = (components...,)
  ntup = NamedTuple{S,typeof(tupcomp)}(tupcomp)
  Mixture(F2, ntup, diversity)
end





  #
