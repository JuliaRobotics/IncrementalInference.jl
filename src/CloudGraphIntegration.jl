

function addCloudVert!(fgl::FactorGraph, exvert::Graphs.ExVertex)
  cv = CloudGraphs.exVertex2CloudVertex(exvert);
  CloudGraphs.add_vertex!(fgl.cg, cv);
  fgl.cgIDs[exvert.index] = cv.neo4jNodeId
  IncrementalInference.addGraphsVert!(fgl, exvert)
end

function makeAddCloudEdge!(fgl::FactorGraph, v1::Graphs.ExVertex, v2::Graphs.ExVertex)
  cv1 = CloudGraphs.get_vertex(fgl.cg, fgl.cgIDs[v1.index], false)
  cv2 = CloudGraphs.get_vertex(fgl.cg, fgl.cgIDs[v2.index], false)
  ce = CloudGraphs.CloudEdge(cv1, cv2, "DEPENDENCE");
  CloudGraphs.add_edge!(fgl.cg, ce);
  IncrementalInference.makeAddEdge!(fgl, v1, v2)
end


function setCloudDataLayerAPI()
  # cgapi = DataLayerAPI(addCloudVert!,            # addvertex
  #                      dlapi.getvertex,          # getvertex
  #                      makeAddCloudEdge!,        # makeaddedge
  #                      dlapi.outneighbors,       # outneighbors
  #                      +, +, +, + )
  dlapi.addvertex! = addCloudVert!
  dlapi.makeaddedge! = makeAddCloudEdge!
  println("Changed internal API calls to use CloudGraphs in appropriate places.")
  nothing
end
