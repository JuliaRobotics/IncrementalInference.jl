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

localapi = DataLayerAPI(addGraphsVert!,          # addvertex
                        getVertNode,             # getvertex
                        makeAddEdge!,            # makeaddedge
                        graphsGetEdge,           # getedge
                        graphsOutNeighbors,      # outneighbors
                        updateFullVertData!,     # updatevertex!
                        +,                       # updateedge!
                        graphsDeleteVertex!,     # deletevertex!
                        +,                       # deleteedge!
                        false )
#


# setCloudDataLayerAPI!
function setdatalayerAPI!(;addvertex::Function = addGraphsVert!,
      getvertex::Function = getVertNode,
      makeaddedge::Function = makeAddEdge!,
      getedge::Function = graphsGetEdge,
      outneighbors::Function = graphsOutNeighbors,
      updatevertex::Function = updateFullVertData!,
      updateedge::Function = +,
      deletevertex::Function = graphsDeleteVertex!,
      deleteedge::Function = +,
      cgEnabled::Bool = false  )

  dlapi.addvertex! = addvertex
  dlapi.getvertex = getvertex
  dlapi.makeaddedge! = makeaddedge
  dlapi.getedge = getedge
  dlapi.outneighbors = outneighbors
  dlapi.updatevertex! = updatevertex
  dlapi.updateedge! = updateedge
  dlapi.deletevertex! = deletevertex
  dlapi.deleteedge! = deleteedge
  dlapi.cgEnabled = cgEnabled

  # dlapi.addvertex! = addCloudVert!
  # dlapi.getvertex = getExVertFromCloud
  # dlapi.makeaddedge! = makeAddCloudEdge!
  # dlapi.getedge = getEdgeFromCloud
  # dlapi.updatevertex! = updateFullCloudVertData!
  # dlapi.outneighbors = getCloudOutNeighbors
  # dlapi.deletevertex! = deleteCloudVertex!
  # dlapi.deleteedge! = deleteCloudEdge!
  # dlapi.cgEnabled = true

  println("Changed internal API calls to use outside calls.")
  nothing
end


function showcurrentdlapi()
  @show dlapi
end





# Remember 3rd party users interact with
# addNode!
# addFactor!
# prepBatchTree!
# inferOverTree!
#
## Visualization functions -- Visualization may well be separated or abstracted out
# writeGraphPdf
# drawHorBeliefsList
