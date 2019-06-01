
"""
    $SIGNATURES

Draw and show the factor graph `<:AbstractDFG` via system graphviz and pdf app.

Notes
- Should not be calling outside programs.
- Need long term solution
- DFG's `toDotFile` a better solution -- view with `xdot` application.
"""
function writeGraphPdf(fgl::G;
                       viewerapp::String="evince",
                       filepath::AS="/tmp/fg.pdf",
                       engine::AS="sfdp",
                       show::Bool=true ) where {G <: AbstractDFG, AS <: AbstractString}
  #
  @warn "writeGraphPdf function might changed, see DFG.toDotFile(dfg) as part of the long term solution."

  fgd = fgl
  @info "Writing factor graph file"
  fext = split(filepath, '.')[end]
  fpwoext = split(filepath, '.')[end-1]
  dotfile = fpwoext*".dot"
  # fid = open(dotfile,"w")
  # write(fid,Graphs.to_dot(fgd.g))
  # close(fid)
  DFG.GraphsJl.toDotFile(fgl, dotfile)
  show ? (@async run(`$(engine) $(dotfile) -T$(fext) -o $(filepath)`)) : nothing

  try
    viewerapp != nothing ? (@async run(`$(viewerapp) $(filepath)`)) : nothing
  catch e
    @warn "not able to show $(filepath) with viewerapp=$(viewerapp). Exception e=$(e)"
  end
  nothing
end


"""
$(SIGNATURES)

Add a node (variable) to a graph. Use this over the other dispatches.

DEPRECATED: use addVarialbe! instead.
"""
function addNode!(fg::FactorGraph,
                  lbl::Symbol,
                  softtype::T;
                  N::Int=100,
                  autoinit::Bool=true,  # does init need to be separate from ready? TODO
                  ready::Int=1,
                  dontmargin::Bool=false,
                  labels::Vector{<:AbstractString}=String[],
                  api::DataLayerAPI=dlapi,
                  uid::Int=-1,
                  smalldata=""  ) where {T <:InferenceVariable}
  #
  @warn "IIF.addNode!(..) is being deprecated, use IIF.addVariable!(..) instead."
  return addVariable!( fg,
                       lbl,
                       softtype,
                       N=N,
                       autoinit=autoinit,  # does init need to be separate from ready? TODO
                       ready=ready,
                       dontmargin=dontmargin,
                       labels=labels,
                       api=api,
                       uid=uid,
                       smalldata=smalldata )
end
function addNode!(fg::FactorGraph,
                  lbl::Symbol,
                  softtype::Type{<:InferenceVariable};
                  N::Int=100,
                  autoinit::Bool=true,
                  ready::Int=1,
                  dontmargin::Bool=false,
                  labels::Vector{<:AbstractString}=String[],
                  api::DataLayerAPI=dlapi,
                  uid::Int=-1,
                  # dims::Int=-1,
                  smalldata=""  )
  #
  @warn "IIF.addNode!(..) is being deprecated, use IIF.addVariable!(..) instead."
  return addVariable!(fg,
                      lbl,
                      softtype,
                      N=N,
                      autoinit=autoinit,
                      ready=ready,
                      dontmargin=dontmargin,
                      labels=labels,
                      api=api,
                      uid=uid,
                      smalldata=smalldata  )
end
