
##==============================================================================
## LEGACY SUPPORT FOR ZMQ IN CAESAR
##==============================================================================

export listSolvekeys

export _evalType

# not sure if and where this is still being used
function _evalType(pt::String)::Type
  @error "_evalType has been deprecated, use DFG serialization methods instead."
  try
    getfield(Main, Symbol(pt))
  catch ex
    io = IOBuffer()
    showerror(io, ex, catch_backtrace())
    err = String(take!(io))
    error("_evalType: Unable to locate factor/distribution type '$pt' in main context (e.g. do a using on all your factor libraries). Please check that this factor type is loaded into main. Stack trace = $err")
  end
end



##==============================================================================
## TODO deprecated  
##==============================================================================


# see DFG #590
@deprecate extractdistribution(x) convert(SamplableBelief, x)


##==============================================================================
## Deprecate at v0.19
##==============================================================================

@deprecate getDwnMsgConsolidated(tree::AbstractBayesTree, edge) getMsgDwnChannel(tree, edge)

# @deprecate putBeliefMessageUp!(tree::AbstractBayesTree, edge, beliefMsg::LikelihoodMessage) putMessageUp!(tree, edge, beliefMsg)
# @deprecate takeBeliefMessageUp!(tree::AbstractBayesTree, edge) takeMessageUp!(tree, edge)
# @deprecate putBeliefMessageDown!(tree::BayesTree, edge, beliefMsg::LikelihoodMessage) putMessageDown!(tree, edge, beliefMsg)
# @deprecate takeBeliefMessageDown!(tree::BayesTree, edge) takeMessageDown!(tree, edge)

@deprecate appendSeparatorToClique(w...;kw...) appendSeparatorToClique!(w...;kw...)

@deprecate TreeClique(i::Int, label::Union{AbstractString, Symbol}) TreeClique(CliqueId(i), BayesTreeNodeData(), Dict{String,Any}())

@deprecate emptyBayesTree() BayesTree()

# Graph.jl does not have an in_edges function for a GenericIncidenceList, so extending here.
# function Graphs.in_edges(vert::V, gr::GenericIncidenceList{V, Edge{V}, Vector{V}}) where {V}
#   inclist = gr.inclist
#   targid = vert.id
#   inlist = Edge{V}[]
#   for edgelist in inclist
#     for ed in edgelist
#       if ed.target.id == targid
#         push!(inlist, ed)
#       end
#     end
#   end
#   return inlist
# end

# getMsgDwnChannel(tree::BayesTree, edge) = tree.messageChannels[edge.index].downMsg

# getMsgUpChannel(tree::BayesTree, edge) = tree.messageChannels[edge.index].upMsg

# function parentCliq(treel::BayesTree, cliq::TreeClique)
#   Graphs.in_neighbors(cliq, treel.bt)
# end
# function parentCliq(treel::BayesTree, frtsym::Symbol)
#   parentCliq(treel,  getClique(treel, frtsym))
# end
# getNumCliqs(tree::BayesTree) = Graphs.num_vertices(tree.bt)

# getEdgesParent(tree::BayesTree, cliq::TreeClique) = Graphs.in_edges(cliq, tree.bt)

# getEdgesChildren(tree::BayesTree, cliq::TreeClique) = Graphs.out_edges(cliq, tree.bt)

# function childCliqs(treel::BayesTree, cliq::TreeClique)
#   childcliqs = Vector{TreeClique}(undef, 0)
#   for cl in Graphs.out_neighbors(cliq, treel.bt)
#     push!(childcliqs, cl)
#   end
#   return childcliqs
# end
# function childCliqs(treel::BayesTree, frtsym::Symbol)
#   childCliqs(treel,  getClique(treel, frtsym))
# end

# function initTreeMessageChannels!(tree::BayesTree)
#   for e = 1:tree.bt.nedges
#     push!(tree.messageChannels, e=>(upMsg=Channel{LikelihoodMessage}(0),downMsg=Channel{LikelihoodMessage}(0)))
#   end
#   return nothing
# end

# isRoot(treel::BayesTree, cliq::TreeClique) = isRoot(tree, cliq.id)
# function isRoot(treel::BayesTree, cliqKey::Int)
#   length(Graphs.in_neighbors(getClique(treel, cliqKey), treel.bt)) == 0
# end

# # TODO Deprecate
# # Graphs.jl BayesTree declarations
# const BTGdict = GenericIncidenceList{TreeClique,Edge{TreeClique},Array{TreeClique,1},Array{Array{Edge{TreeClique},1},1}}

# """
# $(TYPEDEF)

# Data structure for the Bayes (Junction) tree, which is used for inference and constructed from a given `::AbstractDFG`.
# Dev Notes:
# - To be deprecated or `BTGdict` replaced, as Graphs.jl is deprecated.
# """
# mutable struct BayesTree <: AbstractBayesTree
#   bt::BTGdict
#   btid::Int
#   cliques::Dict{Int,TreeClique}
#   frontals::Dict{Symbol,Int}
#   messageChannels::Dict{Int, NamedTuple{(:upMsg, :downMsg),Tuple{Channel{LikelihoodMessage},Channel{LikelihoodMessage}}}}
#   variableOrder::Vector{Symbol}
#   buildTime::Float64
# end

# BayesTree() = BayesTree(Graphs.inclist(TreeClique,is_directed=true),
#                         0,
#                         Dict{Int,TreeClique}(),
#                         Dict{AbstractString, Int}(),
#                         Dict{Int, NamedTuple{(:upMsg, :downMsg),Tuple{Channel{LikelihoodMessage},Channel{LikelihoodMessage}}}}(),
#                         Symbol[],
#                         0.0  )
# #

# getMessageChannels(tree::BayesTree) = tree.messageChannels


# Graphs.make_vertex(g::AbstractGraph{TreeClique}, label::AbstractString) = TreeClique(num_vertices(g) + 1, String(label))
# Graphs.vertex_index(v::TreeClique) = v.id
# Graphs.attributes(v::TreeClique, g::AbstractGraph) = v.attributes

@deprecate initManual!(dfg::AbstractDFG, variable::DFGVariable, ptsArr::BallTreeDensity) initManual!(variable, ptsArr)

@deprecate setCliqueDrawColor!(w...;kw...) setCliqueDrawColor!(w...;kw...)

@deprecate evalFactor2(w...;kw...) evalFactor(w...;kw...)

# TreeBelief field softtype->variableType rename
function Base.getproperty(x::TreeBelief,f::Symbol)
  if f == :softtype
    Base.depwarn("`TreeBelief` field `softtype` is deprecated, use `variableType`", :getproperty)
    f = :variableType
  end
  getfield(x,f)
end

function Base.setproperty!(x::TreeBelief, f::Symbol, val)
  if f == :softtype
    Base.depwarn("`TreeBelief` field `softtype` is deprecated, use `variableType`", :getproperty)
    f = :variableType
  end
  return setfield!(x, f, convert(fieldtype(typeof(x), f), val))
end

# function MetaBayesTree(tree::BayesTree)
#   Base.depwarn("Graphs.jl Bayes Tree is deprecated, this constructor will be removed", :MetaBayesTree)
#   # create graph from Graphs.jl adjacency_matrix
#   mtree = MetaBayesTree(MetaDiGraph{Int, Float64}(MetaGraphs.SimpleDiGraph(Graphs.adjacency_matrix(tree.bt))), tree.btid, tree.frontals, tree.variableOrder, tree.buildTime)

#   #deep copy over properties
#   for v in tree.bt.vertices
#     # set_prop!(mtree.bt, v.id, :label, deepcopy(v.label))
#     set_prop!(mtree.bt, v.id, :clique, deepcopy(v))
#   end

#   return mtree

# end



#
