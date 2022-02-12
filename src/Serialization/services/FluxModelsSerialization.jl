# Serialization functions for Flux models that depend on BSON

# @info "IncrementalInference is adding Flux/BSON serialization functionality."

using Base64

import Base: convert


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


function packDistribution(obj::FluxModelsDistribution)
  #

    # and the specialSampler function -- likely to be deprecated
    # specialSampler = Symbol(obj.specialSampler)
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
    mimeTypeData = "application/octet-stream/bson/base64"
  else
    # store one just model to preserve the type (allows resizing on immutable Ref after deserialize)
    push!(models,_serializeFluxModelBase64(obj.models[1]))
    # at least capture the type of how the data looks for future deserialization
    sdata = string(typeof(obj.data))
    mimeTypeData = "application/text"
  end
  mimeTypeModel = "application/octet-stream/bson/base64"
  
  # and build the JSON-able object
  return PackedFluxModelsDistribution("IncrementalInference.PackedFluxModelsDistribution",
                                      inputDim, 
                                      outputDim, 
                                      mimeTypeModel,
                                      models, 
                                      mimeTypeData,
                                      sdata, 
                                      obj.shuffle[], 
                                      obj.serializeHollow[], 
                                      "IncrementalInference.PackedFluxModelsDistribution" )
  #
end



function unpackDistribution(obj::PackedFluxModelsDistribution)
  #
  obj.serializeHollow && @warn("Deserialization of FluxModelsDistribution.serializationHollow=true is not yet well developed, please open issues at IncrementalInference.jl accordingly.")
  
  # deserialize
  # @assert obj.mimeTypeModel == "application/octet-stream/bson/base64"
  models = _deserializeFluxModelBase64.(obj.models)
  
  # @assert obj.mimeTypeData == "application/octet-stream/bson/base64"
  data = !obj.serializeHollow ? _deserializeFluxDataBase64.(obj.data) : zeros(0)

  return FluxModelsDistribution(models, 
                                (obj.inputDim...,), 
                                data, 
                                (obj.outputDim...,), 
                                shuffle=obj.shuffle, 
                                serializeHollow=obj.serializeHollow  )
end



function Base.convert(::Union{Type{<:PackedSamplableBelief},Type{<:PackedFluxModelsDistribution}}, 
                      obj::FluxModelsDistribution)
  #
  # convert to packed type first
  return packDistribution(obj)
end



function convert( ::Union{Type{<:SamplableBelief},Type{<:FluxModelsDistribution}}, 
                  obj::PackedFluxModelsDistribution )
  #
  return unpackDistribution(obj)
end


#