# API definition file for data layer abstraction.
# The idea is for a one line call to change from internal Graphs.jl to CloudGraphs.jl
# switch for as data layer.

type DataLayerAPI
  addvertex!::Function
  getvertex::Function
  # setupvertgraph!::Function
  # setupfncvertgraph!::Function
  makeaddedge!::Function
  getedge::Function
  outneighbors::Function
  updatevertex!::Function
  updateedge!::Function
  deletevertex!::Function
  deleteedge!::Function
  cgEnabled::Bool
end

# global dlapi

dlapi = DataLayerAPI(addGraphsVert!,          # addvertex
                     getVertNode,             # getvertex
                     makeAddEdge!,            # makeaddedge
                     graphsGetEdge,           # getedge
                     graphsOutNeighbors,      # outneighbors
                     updateFullVertData!,     # updatevertex!
                     +,                       # updateedge!
                     graphsDeleteVertex!,     # deletevertex!
                     +,                       # deleteedge!
                     false )

# Remember 3rd party users interact with
# addNode!
# addFactor!
# prepBatchTree!
# inferOverTree!
#
## Visualization functions -- Visualization may well be separated or abstracted out
# writeGraphPdf
# drawHorBeliefsList
