# deprecated functions

# must set either dims or initval for proper initialization
# Add node to graph, given graph struct, labal, init values,
# std dev [TODO -- generalize], particle size and ready flag for concurrency
function addNode!(fg::FactorGraph,
      lbl::Symbol,
      initval::Array{Float64}=zeros(1,1),
      stdev::Array{Float64}=ones(1,1); # this is bad and should be removed TODO
      N::Int=100,
      ready::Int=1,
      labels::Vector{<: AbstractString}=String[],
      api::DataLayerAPI=dlapi,
      uid::Int=-1,
      dims::Int=-1  )
  #
  warn("this addNode! is deprecated, please use FactorGraph01.jl:addNode!(fg::FactorGraph, lbl::Symbol, softtype::Type{T}) instead.")
  currid = fg.id+1
  if uid==-1
    fg.id=currid
  else
    currid = uid
  end
  dims = dims != -1 ? dims : size(initval,1)

  lblstr = string(lbl)
  vert = ExVertex(currid,lblstr)
  addNewVarVertInGraph!(fg, vert, currid, lbl, ready, nothing)
  # dlapi.setupvertgraph!(fg, vert, currid, lbl) #fg.v[currid]
  dodims = fg.dimID+1
  # TODO -- vert should not loose information here
  setDefaultNodeData!(vert, initval, stdev, dodims, N, dims) #fg.v[currid]

  vnlbls = deepcopy(labels)
  push!(vnlbls, fg.sessionname)
  # addvert!(fg, vert, api=api)
  api.addvertex!(fg, vert, labels=vnlbls) #fg.g ##vertr =

  fg.dimID+=dims # rows indicate dimensions, move to last dimension
  push!(fg.nodeIDs, currid)

  return vert
end

function getVert(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi) where {T <: AbstractString}
  warn("IncrementalInference.getVert{T <: AbstractString}(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi) is deprecated, use lbl::Symbol instead")
  getVert(fgl, Symbol(lbl), api=api)
end

function getVal(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi) where {T <: AbstractString}
  warn("IncrementalInference.getVal{T <: AbstractString}(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi) is deprecated, use lbl::Symbol instead")
  getVal(fgl, Symbol(lbl),api=api)
end
