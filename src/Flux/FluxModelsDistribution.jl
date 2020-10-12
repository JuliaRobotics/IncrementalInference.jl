

@info "IncrementalInference is adding Flux related functionality."



# the factor definitions
export FluxModelsDistribution
# some utilities


# Required packages
using .Flux

using Random, Statistics
# using DistributedFactorGraphs

# import Base: convert
import Random: rand



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




#