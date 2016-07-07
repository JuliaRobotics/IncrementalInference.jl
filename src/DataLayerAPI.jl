# API definition file for data layer abstraction.
# The idea is for a one line call to change from internal Graphs.jl to CloudGraphs.jl
# switch for as data layer.

type DataLayerAPI
  addvertex!::Function
  getvertex::Function
  # setupvertgraph!::Function
  # setupfncvertgraph!::Function
  makeaddedge!::Function
  outneighbors::Function
  updatevertex!::Function
  updateedge!::Function
  deletevertex!::Function
  deleteedge!::Function
end

# global dlapi

dlapi = DataLayerAPI(addGraphsVert!,          # addvertex
                     getVertNode,             # getvertex
                     makeAddEdge!,            # makeaddedge
                     graphsOutNeighbors,      # outneighbors
                     updateFullVertData!,     # updatevertex!
                     +, +, + )

# Remember 3rd party users interact with
# addNode!
# addFactor!
# prepBatchTree!
# inferOverTree!
#
## Visualization functions -- Visualization may well be separated or abstracted out
# writeGraphPdf
# drawHorBeliefsList
