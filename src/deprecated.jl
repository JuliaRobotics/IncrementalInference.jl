# deprecated functions

function getVert(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi) where {T <: AbstractString}
  warn("IncrementalInference.getVert{T <: AbstractString}(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi) is deprecated, use lbl::Symbol instead")
  getVert(fgl, Symbol(lbl), api=api)
end

function getVal(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi) where {T <: AbstractString}
  warn("IncrementalInference.getVal{T <: AbstractString}(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi) is deprecated, use lbl::Symbol instead")
  getVal(fgl, Symbol(lbl),api=api)
end
