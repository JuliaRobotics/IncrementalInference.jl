

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

function getExVertFromCloud(fgl::FactorGraph, lbl::ASCIIString; bigdata::Bool=false)
  getExVertFromCloud(fgl, fgl.IDs[lbl],bigdata=bigdata)
end

function updateFullCloudVertData!(fgl::FactorGraph,
    nv::Graphs.ExVertex)

  # TODO -- this get_vertex seems excessive, but we need the CloudVertex
  neoID =fgl.cgIDs[nv.index]
  # println("updateFullCloudVertData! -- trying to get $(neoID)")
  vert = CloudGraphs.get_vertex(fgl.cg, neoID, false)
  vert.packed = nv.attributes["data"]
  # TODO -- ignoring other properties
  fgl.g.vertices[nv.index].attributes["data"] = nv.attributes["data"]
  CloudGraphs.update_vertex!(fgl.cg, vert)
end

function makeAddCloudEdge!(fgl::FactorGraph, v1::Graphs.ExVertex, v2::Graphs.ExVertex)
  cv1 = CloudGraphs.get_vertex(fgl.cg, fgl.cgIDs[v1.index], false)
  cv2 = CloudGraphs.get_vertex(fgl.cg, fgl.cgIDs[v2.index], false)
  ce = CloudGraphs.CloudEdge(cv1, cv2, "DEPENDENCE");
  retrel = CloudGraphs.add_edge!(fgl.cg, ce);

  # TODO -- keep this edge id in function node data, must refactor
  push!(v2.attributes["data"].edgeIDs, retrel.id) # TODO -- not good way to do this
  updateFullCloudVertData!(fgl, v2)

  IncrementalInference.makeAddEdge!(fgl, v1, v2, saveedgeID=false)
  retrel.id
end

# return list of neighbors as Graphs.ExVertex type
function getCloudOutNeighbors(fgl::FactorGraph, vert::Graphs.ExVertex)
  # @show vert.index
  # @show fgl.cgIDs
  cgid = fgl.cgIDs[vert.index]
  cv = CloudGraphs.get_vertex(fgl.cg, cgid, false)
  neighs = CloudGraphs.get_neighbors(fgl.cg, cv)
  neExV = Graphs.ExVertex[]
  for n in neighs
    push!(neExV,  CloudGraphs.cloudVertex2ExVertex(n))
  end
  return neExV
end

function getEdgeFromCloud(fgl::FactorGraph, id::Int64)
  println("getting id=$(id)")
  CloudGraphs.get_edge(fgl.cg, id)
end

function deleteCloudVertex!(fgl::FactorGraph, vert::Graphs.ExVertex)
  neoID = fgl.cgIDs[vert.index]
  cvr = CloudGraphs.get_vertex(fgl.cg, neoID, false)
  CloudGraphs.delete_vertex!(fgl.cg, cvr)
end

function deleteCloudEdge!(fgl::FactorGraph, edge::CloudEdge)
  CloudGraphs.delete_edge!(fgl.cg, edge)
end

function setCloudDataLayerAPI!()
  # cgapi = DataLayerAPI(addCloudVert!,            # addvertex
  #                      dlapi.getvertex,          # getvertex
  #                      makeAddCloudEdge!,        # makeaddedge
  #                      graphsGetEdge,           # getedge
  #                      dlapi.outneighbors,       # outneighbors
  #                      +, +, +, + )
  dlapi.addvertex! = addCloudVert!
  dlapi.getvertex = getExVertFromCloud
  dlapi.makeaddedge! = makeAddCloudEdge!
  dlapi.getedge = getEdgeFromCloud
  dlapi.updatevertex! = updateFullCloudVertData!
  dlapi.outneighbors = getCloudOutNeighbors
  dlapi.deletevertex! = deleteCloudVertex!
  dlapi.deleteedge! = deleteCloudEdge!
  dlapi.cgEnabled = true

  println("Changed internal API calls to use CloudGraphs in appropriate places.")
  nothing
end






  #
