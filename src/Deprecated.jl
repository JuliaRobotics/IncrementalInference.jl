"""
$(SIGNATURES)

Add a node (variable) to a graph. Use this over the other dispatches.

DEPRECATED: use addVarialbe! instead.
"""
function addNode!(fg::FactorGraph,
                  lbl::Symbol,
                  softtype::T;
                  N::Int=100,
                  autoinit::Bool=true,  # does init need to be separate from ready? TODO
                  ready::Int=1,
                  dontmargin::Bool=false,
                  labels::Vector{<:AbstractString}=String[],
                  api::DataLayerAPI=dlapi,
                  uid::Int=-1,
                  smalldata=""  ) where {T <:InferenceVariable}
  #
  @warn "IIF.addNode!(..) is being deprecated, use IIF.addVariable!(..) instead."
  return addVariable!( fg,
                       lbl,
                       softtype,
                       N=N,
                       autoinit=autoinit,  # does init need to be separate from ready? TODO
                       ready=ready,
                       dontmargin=dontmargin,
                       labels=labels,
                       api=api,
                       uid=uid,
                       smalldata=smalldata )
end
function addNode!(fg::FactorGraph,
                  lbl::Symbol,
                  softtype::Type{<:InferenceVariable};
                  N::Int=100,
                  autoinit::Bool=true,
                  ready::Int=1,
                  dontmargin::Bool=false,
                  labels::Vector{<:AbstractString}=String[],
                  api::DataLayerAPI=dlapi,
                  uid::Int=-1,
                  # dims::Int=-1,
                  smalldata=""  )
  #
  @warn "IIF.addNode!(..) is being deprecated, use IIF.addVariable!(..) instead."
  return addVariable!(fg,
                      lbl,
                      softtype,
                      N=N,
                      autoinit=autoinit,
                      ready=ready,
                      dontmargin=dontmargin,
                      labels=labels,
                      api=api,
                      uid=uid,
                      smalldata=smalldata  )
end
