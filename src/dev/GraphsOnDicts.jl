
using Graphs

# typealias FGGdict Graphs.GenericIncidenceList{Graphs.ExVertex,Graphs.Edge{Graphs.ExVertex},Dict{Int,Graphs.ExVertex},Dict{Int,Array{Graphs.Edge{Graphs.ExVertex},1}}}

# g = Graphs.inclist(Graphs.ExVertex,is_directed=false)
g = incdict(ExVertex,is_directed=false)

v1 = ExVertex(3,"x1")
add_vertex!(g, v1)

v2 = ExVertex(2,"x2")
add_vertex!(g, v2)

e1 = make_edge(g, v1, v2)
add_edge!(g, e1)

vertices(g)

# print neighbors
for n in out_neighbors(v2, g)  @show n.index end

# calculate adjacency_matrix
a, p = adjacency_matrix(g, returnpermutation=true)










#
