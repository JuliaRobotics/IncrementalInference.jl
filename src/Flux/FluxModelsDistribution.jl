

@info "IncrementalInference is adding Flux related functionality."


# the factor definitions
export FluxModelsDistribution
export MixtureFluxModels


# Required packages
using .Flux
using DataStructures: OrderedDict
using Random, Statistics


# import Base: convert
import Random: rand

const _IIFListTypes = Union{<:AbstractVector, <:Tuple, <:NTuple, <:NamedTuple}


function rand(nfb::FluxModelsDistribution,
              N::Int=1 )
  #

  # number of predictors to choose from, and choose random subset
  numModels = length(nfb.models)
  allPreds = 1:numModels |> collect # 1:Npreds |> collect
  # TODO -- compensate when there arent enough prediction models
  if numModels < N
    reps = (N ÷ numModels) + 1
    allPreds = repeat(allPreds, reps  )
    resize!(allPreds, N)
  end
  # samples for the order in which to use models, dont shuffle if N models
  # can suppress shuffle for NN training purposes
  1 < numModels && nfb.shuffle[] ? shuffle!(allPreds) : nothing

  # generate the measurements
  meas = zeros(nfb.outputDim..., N)
  for i in 1:N
    meas[:,i] = (nfb.models[allPreds[i]])(nfb.data)
  end

  return meas
end


FluxModelsDistribution( inDim::NTuple{ID,Int}, 
                        outDim::NTuple{OD,Int}, 
                        models::Vector{P}, 
                        data::D,
                        shuffle::Bool=true,
                        serializeHollow::Bool=false ) where {ID,OD,P,D<:AbstractArray} = FluxModelsDistribution{ID,OD,P,D}(inDim, outDim, models, data, Ref(shuffle), Ref(serializeHollow) )
#

FluxModelsDistribution( models::Vector{P}, 
                        inDim::NTuple{ID,Int}, 
                        data::D,
                        outDim::NTuple{OD,Int};
                        shuffle::Bool=true,
                        serializeHollow::Bool=false ) where {ID,OD,P,D<:AbstractArray} = FluxModelsDistribution{ID,OD,P,D}(inDim, outDim, models, data, Ref(shuffle), Ref(serializeHollow) )
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
function MixtureFluxModels( F_::FunctorInferenceType,
                            nnModels::Vector{P}, 
                            inDim::NTuple{ID,Int}, 
                            data::D,
                            outDim::NTuple{OD,Int},
                            otherComp::_IIFListTypes,
                            diversity::Union{<:AbstractVector, <:NTuple, <:DiscreteNonParametric}; 
                            shuffle::Bool=true,
                            serializeHollow::Bool=false ) where {P,ID,D<:AbstractArray,OD}
  #
  # must preserve order
  allComp = OrderedDict{Symbol, Any}()
  
  # always add the Flux model first
  allComp[:fluxnn] = FluxModelsDistribution(nnModels,
                                            inDim,
                                            data,
                                            outDim,
                                            shuffle=shuffle,
                                            serializeHollow=serializeHollow)
  #
  isNT = otherComp isa NamedTuple
  for idx in 1:length(otherComp)
    nm = isNT ? keys(otherComp)[idx] : Symbol("c$(idx+1)")
    allComp[nm] = otherComp[idx]
  end
  # convert to named tuple
  ntup = (;allComp...)
  
  # construct all the internal objects
  return Mixture(F_, ntup, diversity)
end

MixtureFluxModels(::Type{F}, 
                  w...;
                  kw...) where F <: FunctorInferenceType = MixtureFluxModels(F(LinearAlgebra.I),w...;kw...)





#