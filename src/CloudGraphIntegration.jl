

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
  neoID = fgl.cgIDs[nv.index]
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

  println("Changed internal API calls to use CloudGraphs calls.")
  nothing
end

# register types of interest (Pose2, etc) in CloudGraphs
# you can register new types at any time (Julia is dynamic)
function registerGeneralVariableTypes!(cloudGraph)
  # 1D early development types
  CloudGraphs.registerPackedType!(cloudGraph, VariableNodeData, PackedVariableNodeData, encodingConverter=VNDencoder, decodingConverter=VNDdecoder);
  CloudGraphs.registerPackedType!(cloudGraph, FunctionNodeData{Obsv2}, FunctionNodeData{PackedObsv2}, encodingConverter=FNDencode, decodingConverter=FNDdecode)
  CloudGraphs.registerPackedType!(cloudGraph, FunctionNodeData{Odo}, FunctionNodeData{PackedOdo}, encodingConverter=FNDencode, decodingConverter=FNDdecode)
  CloudGraphs.registerPackedType!(cloudGraph, FunctionNodeData{GenericMarginal}, FunctionNodeData{GenericMarginal}, encodingConverter=FNDencode, decodingConverter=FNDdecode)
  CloudGraphs.registerPackedType!(cloudGraph, FunctionNodeData{Ranged}, FunctionNodeData{Ranged}, encodingConverter=FNDencode, decodingConverter=FNDdecode)
  # Pose2
  CloudGraphs.registerPackedType!(cloudGraph, FunctionNodeData{PriorPose2}, FunctionNodeData{PackedPriorPose2}, encodingConverter=FNDencode, decodingConverter=FNDdecode)
  CloudGraphs.registerPackedType!(cloudGraph, FunctionNodeData{Pose2Pose2}, FunctionNodeData{PackedPose2Pose2}, encodingConverter=FNDencode, decodingConverter=FNDdecode)
  CloudGraphs.registerPackedType!(cloudGraph, FunctionNodeData{Pose2DPoint2DBearingRange}, FunctionNodeData{PackedPose2DPoint2DBearingRange}, encodingConverter=FNDencode, decodingConverter=FNDdecode)
  # TODO -- Pose3 stuff
  nothing
end


# function should not be necessary, but fixes a minor bug following elimination algorithm
function removeGenericMarginals!(conn)
  loadtx = transaction(conn)
  query = "match (n) where n.packedType = 'IncrementalInference.FunctionNodeData{IncrementalInference.GenericMarginal}' delete n"
  cph = loadtx(query, submit=true)
  loadresult = commit(loadtx)
  nothing
end


function getAllExVertexNeoIDs(conn)

  loadtx = transaction(conn)
  query = "match (n) where n.ready=1 and n.backendset=1 return n"
  cph = loadtx(query, submit=true)
  ret = Array{Tuple{Int64,Int64},1}()

  for data in cph.results[1]["data"]
    metadata = data["meta"][1]
    rowdata = data["row"][1]
    push!(ret, (rowdata["exVertexId"],metadata["id"])  )
  end
  return ret
end

# function getDBAdjMatrix()
#
# end

function copyAllEdges!(fgl::FactorGraph, cverts::Dict{Int64, CloudVertex}, IDs::Array{Tuple{Int64,Int64},1})
  # do entire graph, one node at a time
  for ids in IDs
    # look at neighbors of this node
    for nei in CloudGraphs.get_neighbors(fgl.cg, cverts[ids[2]])
      # want to ignore if the edge was previously added from the other side, comparing to the out neighbors in the Graphs structure

      if nei.properties["ready"]==1
          alreadythere = false
          v2 = fgl.g.vertices[nei.exVertexId]
          for graphsnei in Graphs.out_neighbors(v2, fgl.g)
            if graphsnei.index == nei.exVertexId
              alreadythere = true
              break;
            end
          end
          if !alreadythere
            # add the edge to graph
            v1 = fgl.g.vertices[ids[1]]
            makeAddEdge!(fgl, v1, v2, saveedgeID=false)
          end
      end


    end
  end
  nothing
end

function copyAllNodes!(fgl::FactorGraph, cverts::Dict{Int64, CloudVertex}, IDs::Array{Tuple{Int64,Int64},1}, conn)
  for ids in IDs
    cvert = CloudGraphs.get_vertex(fgl.cg, ids[2], false)
    cverts[ids[2]] = cvert
    exvert = cloudVertex2ExVertex(cvert)
    Graphs.add_vertex!(fgl.g, exvert)
    fgl.id < exvert.index ? fgl.id = exvert.index : nothing
    fgl.cgIDs[ids[1]] = ids[2]
    if typeof(exvert.attributes["data"]) == VariableNodeData  # variable node
      fgl.IDs[exvert.label] = ids[1]
      push!(fgl.nodeIDs, ids[1])
    else # function node
      fgl.fIDs[exvert.label] = ids[1]
      push!(fgl.factorIDs, ids[1])
    end
  end
  nothing
end

function fullLocalGraphCopy!(fgl::FactorGraph, conn)

  IDs = getAllExVertexNeoIDs(conn)
  if length(IDs) > 1
    cverts = Dict{Int64, CloudVertex}()
    unsorted = Int64[]
    # TODO ensure this is row is sorted
    for ids in IDs push!(unsorted, ids[1]) end
    testlist = deepcopy(unsorted)
    if testlist != sort(unsorted)
      # TODO -- maybe not required, but being safe for now
      error("Must be sorted list for elimination...")
    end

    # get and add all the nodes
    copyAllNodes!(fgl, cverts, IDs, conn)

    # get and insert all edges
    copyAllEdges!(fgl, cverts, IDs)
    return true
  else
    print(".")
    return false
  end
end

function setDBAllReady!(conn)
  loadtx = transaction(conn)
  query = "match (n) where n.ready=0 set n.ready=1"
  cph = loadtx(query, submit=true)
  loadresult = commit(loadtx)
  nothing
end

function setBackendWorkingSet!(conn)
  loadtx = transaction(conn)
  query = "match (n) set n.backendset=1"
  cph = loadtx(query, submit=true)
  loadresult = commit(loadtx)
  nothing
end








  #
