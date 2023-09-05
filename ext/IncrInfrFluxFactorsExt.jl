module IncrInfrFluxFactorsExt


@info "IncrementalInference is loading extension functionality related to Flux.jl"

# Required packages
using Flux
using DataStructures: OrderedDict
using LinearAlgebra
using Base64
using Manifolds
using DocStringExtensions
using BSON

import Base: convert

# import Base: convert
using Random, Statistics
import Random: rand

using IncrementalInference
import IncrementalInference: samplePoint, sampleTangent, MixtureFluxModels, getSample

# the factor definitions
# export FluxModelsDistribution
export MixtureFluxModels

const _IIFListTypes = Union{<:AbstractVector, <:Tuple, <:NTuple, <:NamedTuple}

function Random.rand(nfb::FluxModelsDistribution, N::Integer = 1)
  #

  # number of predictors to choose from, and choose random subset
  numModels = length(nfb.models)
  allPreds = 1:numModels |> collect
  # TODO -- compensate when there arent enough prediction models
  if !(N isa Nothing) && numModels < N
    reps = (N ÷ numModels) + 1
    allPreds = repeat(allPreds, reps)
    resize!(allPreds, N)
  end
  # samples for the order in which to use models, dont shuffle if N models
  # can suppress shuffle for NN training purposes
  selPred = 1 < numModels && nfb.shuffle[] ? rand(allPreds, N) : view(allPreds, 1:N)

  # dev function, TODO simplify to direct call 
  _sample() = map(pred -> (nfb.models[pred])(nfb.data), selPred)

  return _sample()
  # return [_sample() for _ in 1:N]
end

sampleTangent(M::AbstractManifold, fmd::FluxModelsDistribution, p = 0) = rand(fmd, 1)[1]
samplePoint(M::AbstractManifold, fmd::FluxModelsDistribution, p = 0) = rand(fmd, 1)[1]
function samplePoint(M::AbstractDecoratorManifold, fmd::FluxModelsDistribution, p = 0)
  return rand(fmd, 1)[1]
end

function FluxModelsDistribution(
  inDim::NTuple{ID, Int},
  outDim::NTuple{OD, Int},
  models::Vector{P},
  data::D,
  shuffle::Bool = true,
  serializeHollow::Bool = false,
) where {ID, OD, P, D <: AbstractArray}
  return FluxModelsDistribution{ID, OD, P, D}(
    inDim,
    outDim,
    models,
    data,
    Ref(shuffle),
    Ref(serializeHollow),
  )
end
#

function FluxModelsDistribution(
  models::Vector{P},
  inDim::NTuple{ID, Int},
  data::D,
  outDim::NTuple{OD, Int};
  shuffle::Bool = true,
  serializeHollow::Bool = false,
) where {ID, OD, P, D <: AbstractArray}
  return FluxModelsDistribution{ID, OD, P, D}(
    inDim,
    outDim,
    models,
    data,
    Ref(shuffle),
    Ref(serializeHollow),
  )
end
#

"""
    $SIGNATURES

Helper function to construct `MixtureFluxModels` containing a `NamedTuple`, resulting in a 
`::Mixture` such that `(fluxnn=FluxNNModels, c1=>MvNormal, c2=>Uniform...)` and order sensitive 
`diversity=[0.7;0.2;0.1]`.  The result is the mixture heavily favors `.fluxnn` and names 
`c1` and `c2` for two other components were auto generated.

Notes
- The user can specify own component names if desired (see example).
- `shuffle` is passed through to internal `FluxModelsDistribution` to command shuffling of NN models.
  - `shuffle` does not influence selection of components in the mixture.

Example:

```julia
# some made up data
data = randn(10)
# Flux models
models = [Flux.Chain(softmax, Dense(10,5,σ), Dense(5,1, tanh)) for i in 1:20]
# mixture with user defined names (optional) -- could also just pass Vector or Tuple of components
mix = MixtureFluxModels(PriorSphere1, models, (10,), data, (1,), 
                        (naiveNorm=Normal(),naiveUnif=Uniform()),
                        [0.7; 0.2; 0.1],
                        shuffle=false )
#

# test by add to simple graph
fg = initfg()
addVariable!(fg, :testmix, Sphere1)
addFactor!(fg, [:testmix;], mix)

# look at proposal distribution from the only factor on :testmix
_,pts,__, = localProduct(fg, :testmix)
```

Related

Mixture, FluxModelsDistribution
"""
function MixtureFluxModels(
  F_::AbstractFactor,
  nnModels::Vector{P},
  inDim::NTuple{ID, Int},
  data::D,
  outDim::NTuple{OD, Int},
  otherComp::_IIFListTypes,
  diversity::Union{<:AbstractVector, <:NTuple, <:DiscreteNonParametric};
  shuffle::Bool = true,
  serializeHollow::Bool = false,
) where {P, ID, D <: AbstractArray, OD}
  #
  # must preserve order
  allComp = OrderedDict{Symbol, Any}()

  # always add the Flux model first
  allComp[:fluxnn] = FluxModelsDistribution(
    nnModels,
    inDim,
    data,
    outDim;
    shuffle = shuffle,
    serializeHollow = serializeHollow,
  )
  #
  isNT = otherComp isa NamedTuple
  for idx = 1:length(otherComp)
    nm = isNT ? keys(otherComp)[idx] : Symbol("c$(idx+1)")
    allComp[nm] = otherComp[idx]
  end
  # convert to named tuple
  ntup = (; allComp...)

  # construct all the internal objects
  return Mixture(F_, ntup, diversity)
end

function MixtureFluxModels(::Type{F}, w...; kw...) where {F <: AbstractFactor}
  return MixtureFluxModels(F(LinearAlgebra.I), w...; kw...)
end

#

include("FluxModelsSerialization.jl")


end # module