# API definition file for data layer abstraction.
# The idea is for a one line call to change from internal Graphs.jl to CloudGraphs.jl
# switch for as data layer.

type DataLayerAPI
  addvertex!::Function
  getvertex::Function
  setupvertgraph!::Function
  setupfncvertgraph!::Function
  makeedge::Function
  addedge!::Function
  outneighbors::Function
  updatevertex!::Function
  updateedge!::Function
  deletevertex!::Function
  deleteedge!::Function
end

dlapi = DataLayerAPI(Graphs.add_vertex!,   # addvertex
                    getVertNode,   # getvertex
                    addNewVarVertInGraph!,   # setupvertgraph
                    addNewFncVertInGraph!,   # setupfncvertgraph
                    Graphs.make_edge,   # makeedge
                    Graphs.add_edge!,   # addedge
                    Graphs.out_neighbors,   # outneighbors
                    +, +, +, + )

function setDataLayerAPI(dl::DataLayerAPI)
  IncrementalInference.dlapi = dl
  nothing
end
