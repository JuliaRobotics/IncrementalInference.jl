using IncrementalInference, CloudGraphs, Neo4j

# switch IncrementalInference to use CloudGraphs (Neo4j) data layer

# connect to the server, CloudGraph stuff
configuration = CloudGraphs.CloudGraphConfiguration("localhost", 7474, "", "", "localhost", 27017, false, "", "");

cloudGraph = connect(configuration);

# register types of interest in CloudGraphs
registerGeneralVariableTypes!(cloudGraph)

IncrementalInference.setCloudDataLayerAPI!()

# Connect to database
conn = cloudGraph.neo4j.connection

# function should not be necessary, but fixes a minor bug following elimination algorithm
removeGenericMarginals!(conn)

# TODO -- MAKE INCREMENAL in graph, SUBGRAPHS work in progress!!!!!
while true
  setBackendWorkingSet!(conn)
  # this is being replaced by cloudGraph, added here for development period
  fg = emptyFactorGraph()
  fg.cg = cloudGraph

  # okay now do the solve
  if fullLocalGraphCopy!(fg, conn)
    tree = wipeBuildNewTree!(fg)
    removeGenericMarginals!(conn)
    inferOverTree!(fg, tree)
  else
    sleep(0.2)
  end
end




  #
