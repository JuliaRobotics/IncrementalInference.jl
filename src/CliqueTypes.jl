# Clique types


## Cliques

"""
    $(TYPEDEF)
Structure to store clique data
DEV NOTES: To replace TreeClique completely
    $(FIELDS)
"""
mutable struct TreeClique
  index::Int # see issue #540
  label::Symbol #NOTE this is currently a label such as clique 1, # The drawing label is saved in attributes, JT I'm not sure of the current use
  data::Any#BayesTreeNodeData #FIXME There is circular type usage in TreeClique, BayesTreeNodeData, CliqStateMachineContainer https://github.com/JuliaLang/julia/issues/269
  attributes::Dict{String, Any} #The drawing attributes
  #solveInProgress #on a clique level a "solve in progress" might be very handy
end

TreeClique(i::Int, label::Symbol) = TreeClique(i, label, BayesTreeNodeData(), Dict{String,Any}())
TreeClique(i::Int, label::AbstractString) = TreeClique(i, Symbol(label))

Graphs.make_vertex(g::AbstractGraph{TreeClique}, label::AbstractString) = TreeClique(num_vertices(g) + 1, String(label))
Graphs.vertex_index(v::TreeClique) = v.index
Graphs.attributes(v::TreeClique, g::AbstractGraph) = v.attributes

#TODO the label field and label atribute is a bit confusing with accessors.
DFG.getLabel(cliq::TreeClique) = cliq.attributes["label"]
function setLabel!(cliq::TreeClique, lbl::String)
  cliq.attributes["label"] = lbl
  lbl
end


function compare(c1::TreeClique,
                 c2::TreeClique )
  #
  TP = true
  TP = TP && c1.index == c2.index
  TP = TP && c1.label == c2.label
  # data
  @warn "skipping ::TreeClique compare of data"
  # TP = TP && compare(c1.data, c2.data)

  # attributes
  @warn "only comparing keys of TreeClique attributes"
  TP = TP && collect(keys(c1.attributes)) == collect(keys(c2.attributes))

  return TP
end

## end Cliques
