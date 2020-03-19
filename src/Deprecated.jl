
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


# function writeGraphPdf(fgl::G;
#                        viewerapp::String="evince",
#                        filepath::AS="/tmp/fg.pdf",
#                        engine::AS="neato", #sfdp
#                        show::Bool=true ) where {G <: AbstractDFG, AS <: AbstractString}
#   #
#   @warn "writeGraphPdf is function changing to drawGraph, see DFG.toDotFile(dfg) as part of the long term solution."
#
#   fgd = fgl
#   @info "Writing factor graph file"
#   fext = split(filepath, '.')[end]
#   fpwoext = filepath[1:(end-length(fext)-1)] # split(filepath, '.')[end-1]
#   dotfile = fpwoext*".dot"
#
#   # create the dot file
#   DFG.toDotFile(fgl, dotfile)
#
#   try
#     run(`$(engine) $(dotfile) -T$(fext) -o $(filepath)`)
#     show ? (@async run(`$(viewerapp) $(filepath)`)) : nothing
#   catch e
#     @warn "not able to show $(filepath) with viewerapp=$(viewerapp). Exception e=$(e)"
#   end
#   nothing
# end




##==============================================================================
## Delete in v0.10.x if possible, but definitely by v0.11
##==============================================================================

@deprecate setData!(v::TreeClique, data) setCliqueData!(v,data)


##==============================================================================
## Delete in v0.11
##==============================================================================

"""
    $SIGNATURES

writeGraphPdf deprecated, use drawGraph instead
"""
function writeGraphPdf(fgl::AbstractDFG;
                       viewerapp::AbstractString="evince",
                       filepath::AbstractString="/tmp/fg.pdf",
                       engine::AbstractString="neato",
                       show::Bool=true )
  #
  @warn "writeGraphPdf deprecated, use drawGraph instead"
  drawGraph(fgl, viewerapp=viewerapp, filepath=filepath, engine=engine, show=show )
end

@deprecate manualinit!(dfg::AbstractDFG, vert::DFGVariable, pX::BallTreeDensity) initManual!(dfg::AbstractDFG, vert::DFGVariable, pX::BallTreeDensity)
@deprecate manualinit!(dfg::AbstractDFG, sym::Symbol, pX::BallTreeDensity) initManual!(dfg::AbstractDFG, sym::Symbol, pX::BallTreeDensity) false
@deprecate manualinit!(dfg::AbstractDFG, sym::Symbol, usefcts::Vector{Symbol}) initManual!(dfg::AbstractDFG, sym::Symbol, usefcts::Vector{Symbol}) false
@deprecate manualinit!(dfg::AbstractDFG, sym::Symbol, pts::Array{Float64,2}) initManual!(dfg::AbstractDFG, sym::Symbol, pts::Array{Float64,2}) false


#
