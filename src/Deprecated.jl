
### DONT DELETE YET -- see more likely list below


# function getShortestPathNeighbors(fgl::FactorGraph;
#     from::TreeClique=nothing,
#     to::TreeClique=nothing,
#     neighbors::Int=0 )
#
#   edgelist = shortest_path(fgl.g, ones(num_edges(fgl.g)), from, to)
#   vertdict = Dict{Int,TreeClique}()
#   edgedict = edgelist2edgedict(edgelist)
#   expandVertexList!(fgl, edgedict, vertdict) # grow verts
#   for i in 1:neighbors
#     expandEdgeListNeigh!(fgl, vertdict, edgedict) # grow edges
#     expandVertexList!(fgl, edgedict, vertdict) # grow verts
#   end
#   return vertdict
# end

# function subgraphShortestPath(fgl::FactorGraph;
#                               from::TreeClique=nothing,
#                               to::TreeClique=nothing,
#                               neighbors::Int=0  )
#   #
#   vertdict = getShortestPathNeighbors(fgl, from=from, to=to, neighbors=neighbors)
#   return genSubgraph(fgl, vertdict)
# end



"""
    $SIGNATURES

Draw and show the factor graph `<:AbstractDFG` via system graphviz and pdf app.

Notes
- Should not be calling outside programs.
- Need long term solution
- DFG's `toDotFile` a better solution -- view with `xdot` application.
- also try `engine={"sfdp","fdp","dot","twopi","circo","neato"}`

Future:
- Might be kept with different call strategy since this function so VERY useful!
- Major issue that this function calls an external program such as "evince", which should be
   under user control only.
- Maybe solution is
- ```toDot(fg,file=...); @async run(`xdot file.dot`)```, or
  - ```toDot(fg,file=...); exportPdf(...); @async run(`evince ...pdf`)```.
"""
function writeGraphPdf(fgl::G;
                       viewerapp::String="evince",
                       filepath::AS="/tmp/fg.pdf",
                       engine::AS="neato", #sfdp
                       show::Bool=true ) where {G <: AbstractDFG, AS <: AbstractString}
  #
  @warn "writeGraphPdf is function changing to drawGraph, see DFG.toDotFile(dfg) as part of the long term solution."

  fgd = fgl
  @info "Writing factor graph file"
  fext = split(filepath, '.')[end]
  fpwoext = filepath[1:(end-length(fext)-1)] # split(filepath, '.')[end-1]
  dotfile = fpwoext*".dot"

  # create the dot file
  DFG.toDotFile(fgl, dotfile)

  try
    run(`$(engine) $(dotfile) -T$(fext) -o $(filepath)`)
    show ? (@async run(`$(viewerapp) $(filepath)`)) : nothing
  catch e
    @warn "not able to show $(filepath) with viewerapp=$(viewerapp). Exception e=$(e)"
  end
  nothing
end



##==============================================================================
## Delete in v0.10.x if possible, but definitely by v0.11
##==============================================================================



# For Bayes tree
function setData!(v::TreeClique, data)
  @warn "IIF.setData! should not be used, use IIF.setCliqueData!"
  # this is a memory gulp without replacement, old attr["data"] object is left to gc
  # v.attributes["data"] = data
  # error("dont IIF.setData! for TreeClique")
  v.data = data
  nothing
end




##==============================================================================
## Delete in v0.11
##==============================================================================


@deprecate manualinit!(dfg::AbstractDFG, vert::DFGVariable, pX::BallTreeDensity) initManual!(dfg::AbstractDFG, vert::DFGVariable, pX::BallTreeDensity)
@deprecate manualinit!(dfg::AbstractDFG, sym::Symbol, pX::BallTreeDensity) initManual!(dfg::AbstractDFG, sym::Symbol, pX::BallTreeDensity) false
@deprecate manualinit!(dfg::AbstractDFG, sym::Symbol, usefcts::Vector{Symbol}) initManual!(dfg::AbstractDFG, sym::Symbol, usefcts::Vector{Symbol}) false
@deprecate manualinit!(dfg::AbstractDFG, sym::Symbol, pts::Array{Float64,2}) initManual!(dfg::AbstractDFG, sym::Symbol, pts::Array{Float64,2}) false


#
