# deprecated functions

function getVert{T <: AbstractString}(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi)
  warn("IncrementalInference.getVert{T <: AbstractString}(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi) is deprecated, use lbl::Symbol instead")
  getVert(fgl, Symbol(lbl), api=api)
end

function getVal{T <: AbstractString}(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi)
  warn("IncrementalInference.getVal{T <: AbstractString}(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi) is deprecated, use lbl::Symbol instead")
  getVal(fgl, Symbol(lbl),api=api)
end
