

function addCloudVert!(fgl::FactorGraph, exvert::Graphs.ExVertex)
  cv = CloudGraphs.exVertex2CloudVertex(exvert);
  CloudGraphs.add_vertex!(fgl.cg, cv);
  fgl.cgIDs[exvert.index] = cv.neo4jNodeId
  IncrementalInference.addGraphsVert!(fgl, exvert)
end

# Return Graphs.ExVertex type containing data according to id
function getExVertFromCloud(fgl::FactorGraph, fgid::Int64; bigdata::Bool=false)
  neoID = fgl.cgIDs[fgid]
  cvr = CloudGraphs.get_vertex(fgl.cg, neoID, false)
  CloudGraphs.cloudVertex2ExVertex(cvr)
end

function updateFullCloudVertData!(fgl::FactorGraph,
    nv::Graphs.ExVertex)

  # TODO -- this get_vertex seems excessive, but we need the CloudVertex
  @show nv.index
  neoID =fgl.cgIDs[nv.index]
  println("updateFullCloudVertData! -- trying to get $(neoID)")
  vert = CloudGraphs.get_vertex(fgl.cg, neoID, false)
  vert.packed = nv.attributes["data"]
  # TODO -- ignoring other properties
  CloudGraphs.update_vertex!(fgl.cg, vert)
end

function makeAddCloudEdge!(fgl::FactorGraph, v1::Graphs.ExVertex, v2::Graphs.ExVertex)
  cv1 = CloudGraphs.get_vertex(fgl.cg, fgl.cgIDs[v1.index], false)
  cv2 = CloudGraphs.get_vertex(fgl.cg, fgl.cgIDs[v2.index], false)
  ce = CloudGraphs.CloudEdge(cv1, cv2, "DEPENDENCE");
  CloudGraphs.add_edge!(fgl.cg, ce);
  IncrementalInference.makeAddEdge!(fgl, v1, v2)
end

# return list of neighbors as Graphs.ExVertex type
function getCloudOutNeighbors(fgl::FactorGraph, vert::Graphs.ExVertex)
  @show vert.index
  @show fgl.cgIDs
  cgid = fgl.cgIDs[vert.index]
  cv = CloudGraphs.get_vertex(fgl.cg, cgid, false)
  neighs = CloudGraphs.get_neighbors(fgl.cg, cv)
  neExV = Graphs.ExVertex[]
  for n in neighs
    push!(neExV,  CloudGraphs.cloudVertex2ExVertex(n))
  end
  return neExV
end

function setCloudDataLayerAPI()
  # cgapi = DataLayerAPI(addCloudVert!,            # addvertex
  #                      dlapi.getvertex,          # getvertex
  #                      makeAddCloudEdge!,        # makeaddedge
  #                      dlapi.outneighbors,       # outneighbors
  #                      +, +, +, + )
  dlapi.addvertex! = addCloudVert!
  dlapi.getvertex = getExVertFromCloud
  dlapi.makeaddedge! = makeAddCloudEdge!
  dlapi.updatevertex! = updateFullCloudVertData!
  dlapi.outneighbors = getCloudOutNeighbors
  println("Changed internal API calls to use CloudGraphs in appropriate places.")
  nothing
end
