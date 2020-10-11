# Serialization functions for Flux models that depend on BSON

@info "IncrementalInference is adding Flux/BSON serialization functionality."

export PackedFluxModelsDistribution

using .BSON
using Base64

import Base: convert


mutable struct PackedFluxModelsDistribution <: IIF.PackedSamplableBelief
  # shape of the input data
  inputDim::Vector{Int}
  # shape of the output data
  outputDim::Vector{Int}
  # actual Flux models (Base64 encoded binary)
  mimeTypeModel::String
  models::Vector{String}
  # the data used for prediction, must be <: AbstractArray
  mimeTypeData::String
  data::String
  # shuffle model predictions relative to particle index at each sampling
  shuffle::Bool
  # false for default serialization with model info, set true for separate storage of models 
  serializeHollow::Bool
  # TODO remove requirement and standardize sampler API
  specialSampler::Symbol
  # field name usage to direct the IIF serialization towards JSON method
  PackedSamplableTypeJSON::String
end

function _serializeFluxModelBase64(model::Flux.Chain)
  io = IOBuffer()
  iob64 = Base64EncodePipe(io)
  BSON.@save iob64 model
  close(iob64)
  return String(take!(io))
end

function _deserializeFluxModelBase64(smodel::AbstractString)
  iob64 = PipeBuffer(base64decode(smodel))
  BSON.@load iob64 model
  close(iob64)
  return model
end

function _serializeFluxDataBase64(data::AbstractArray)
  io = IOBuffer()
  iob64 = Base64EncodePipe(io)
  BSON.@save iob64 data
  close(iob64)
  return String(take!(io))
end

function _deserializeFluxDataBase64(sdata::AbstractString)
  iob64 = PipeBuffer(base64decode(sdata))
  BSON.@load iob64 data
  close(iob64)
  return data
end


function convert( ::Union{Type{<:PackedSamplableBelief},Type{<:PackedFluxModelsDistribution}}, 
                  obj::FluxModelsDistribution)
  #

  # and the specialSampler function -- likely to be deprecated
  specialSampler = Symbol(obj.specialSampler)
  # fields to persist
  inputDim = collect(obj.inputDim)
  outputDim = collect(obj.outputDim)
  models = Vector{String}()
  # store all models as Base64 Strings (using BSON)
  if !obj.serializeHollow[]
    resize!(models, length(obj.models))
    # serialize the Vector of Flux models (each one individually)
    models .= _serializeFluxModelBase64.(obj.models)
    # also store data as Base64 String, using BSON
    sdata = _serializeFluxDataBase64(obj.data)
    mimeTypeData = "application/bson/octet-stream/base64"
  else
    # store one just model to preserve the type (allows resizing on immutable Ref after deserialize)
    push!(models,_serializeFluxModelBase64(obj.models[1]))
    # at least capture the type of how the data looks for future deserialization
    sdata = string(typeof(obj.data))
    mimeTypeData = "application/text"
  end
  mimeTypeModel = "application/bson/octet-stream/base64"
  
  # and build the JSON-able object
  packed = PackedFluxModelsDistribution(inputDim, 
                                        outputDim, 
                                        mimeTypeModel,
                                        models, 
                                        mimeTypeData,
                                        sdata, 
                                        obj.shuffle[], 
                                        obj.serializeHollow[], 
                                        specialSampler,
                                        "IncrementalInference.PackedFluxModelsDistribution" )
  #
  return JSON2.write(packed)
end



function convert( ::Union{Type{<:SamplableBelief},Type{FluxModelsDistribution}}, 
                  obj::PackedFluxModelsDistribution)
  #

  obj.serializeHollow && @warn("Deserialization of FluxModelsDistribution.serializationHollow=true is not yet well developed, please open issues at IncrementalInference.jl accordingly.")

  # specialSampler likely to be deprecated
  # specialSampler = getfield(Main, obj.specialSampler)
  
  # deserialize
  # @assert obj.mimeTypeModel == "application/bson/octet-stream/base64"
  models = _deserializeFluxModelBase64.(obj.models)
  
  # @assert obj.mimeTypeData == "application/bson/octet-stream/base64"
  data = !obj.serializeHollow ? _deserializeFluxDataBase64.(obj.data) : zeros(0)

  return FluxModelsDistribution((obj.inputDim...,), 
                                (obj.outputDim...,), 
                                models, 
                                data, 
                                obj.shuffle, 
                                obj.serializeHollow  )
end


#