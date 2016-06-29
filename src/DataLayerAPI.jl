# API definition file for data layer abstraction.
# The idea is for a one line call to change from internal Graphs.jl to CloudGraphs.jl
# switch for as data layer.

type DataLayerAPI
  addvertex!::Function
  getvertex::Function
  setupvertgraph!::Function
  makeedge::Function
  addedge!::Function
  outneighbors::Function
  updatevertex!::Function
  updateedge!::Function
  deletevertex!::Function
  deleteedge!::Function
end

dlapi = DataLayerAPI(Graphs.add_vertex!,
                    getVarNode,
                    addNewVertInGraph!,
                    Graphs.make_edge,
                    Graphs.add_edge!,
                    Graphs.out_neighbors,
                    +, +, +, + )

function setDataLayerAPI(dl::DataLayerAPI)
  IncrementalInference.dlapi = dl
  nothing
end
