
"""
    $SIGNATURES

Construct (new) subgraph and draw the subgraph associated with clique `frontalSym::Symbol`.

Notes
- See `drawGraphCliq`/`writeGraphPdf` for details on keyword options.

Related

drawGraphCliq, spyCliqMat, drawTree, buildCliqSubgraphUp, buildSubgraphFromLabels!
"""
function drawCliqSubgraphUpMocking(
  fgl::G,
  treel::AbstractBayesTree,
  frontalSym::Symbol;
  show::Bool = true,
  filepath::String = "/tmp/caesar/random/cliq_sfg.dot",
  engine::AS1 = "sfdp",
  viewerapp::AS2 = "xdot",
) where {G <: AbstractDFG, AS1 <: AbstractString, AS2 <: AbstractString}
  #
  sfg = buildCliqSubgraphUp(fgl, treel, frontalSym)
  drawGraph(sfg; show = show, viewerapp = viewerapp, engine = engine, filepath = filepath)
  return nothing
end

"""
    $SIGNATURES

Draw and show the factor graph `<:AbstractDFG` via system graphviz and xdot app.

Notes
- Requires system install on Linux of `sudo apt-get install xdot`
- Should not be calling outside programs.
- Need long term solution
- DFG's `DFG.toDotFile` a better solution -- view with `xdot` application.
- also try `engine={"sfdp","fdp","dot","twopi","circo","neato"}`

Notes:
- Calls external system application `xdot` to read the `.dot` file format
  - ```DFG.toDot(fg,file=...); @async run(`xdot file.dot`)```

Related

drawGraphCliq, [`drawTree`](@ref), printCliqSummary, spyCliqMat
"""
function drawGraph(
  fgl::AbstractDFG;
  viewerapp::AbstractString = "xdot",
  filepath::AbstractString = "/tmp/caesar/random/fg.dot",
  engine::AbstractString = "neato", #sfdp
  show::Bool = true,
)
  #
  mkpath(dirname(filepath))
  #   mkpath(joinpath( "/", (split(filepath, '/')[1:(end-1)])...) )  

  @debug "Writing factor graph file"
  fext = split(filepath, '.')[end]
  fpwoext = filepath[1:(end - length(fext) - 1)] # split(filepath, '.')[end-1]
  dotfile = fpwoext * ".dot"

  # create the dot file
  DFG.toDotFile(fgl, dotfile)

  try
    # run(`$(engine) $(dotfile) -T$(fext) -o $(filepath)`)
    show ? (@async run(`$(viewerapp) $(dotfile)`)) : nothing
  catch e
    @warn "not able to show $(filepath) with viewerapp=$(viewerapp). Exception e=$(e)"
  end
  return nothing
end

"""
    $SIGNATURES

Draw the factor graph from a clique state machine history at a particular step as pdf and show.

Related

drawCliqSubgraphUpMocking, drawGraph, drawTree
"""
function drawGraphCliq(
  hists::Dict{Int, <:Tuple},
  step::Int,
  tree::AbstractBayesTree,
  frontal::Symbol;
  show::Bool = true,
)
  #
  cid = getId(getClique(tree, frontal))
  cfg = hists[cid][step][4].cliqSubFg
  return drawGraph(cfg; show = show)
end

"""
    $SIGNATURES

Print basic statistics about a clique variables and factors.

Related

printCliqHistorySummary
"""
function printCliqSummary(
  dfg::G,
  cliq::TreeClique,
  logger = ConsoleLogger(),
) where {G <: AbstractDFG}
  #
  frtl = getCliqFrontalVarIds(cliq)
  seps = getCliqSeparatorVarIds(cliq)
  fcts = getCliqFactorIdsAll(cliq)

  isinit = map(x -> isInitialized(dfg, x), [frtl; seps])
  # infdim = map(x->getVariableInferredDim(dfg, x), [frtl;seps])

  with_logger(logger) do
    @info "Clique $(getId(cliq)) summary:"
    @info "  num frontals:    $(length(frtl))"
    @info "  num separators:  $(length(seps))"
    @info "  num factors:     $(length(fcts))"
    @info "  num initialized: $(sum(isinit)) of $(length(isinit))"
    @info ""
    @info "  frontals:  $(frtl)"
    @info "  separator: $(seps)"
    @info "  factors:   $(fcts)"
    @info "  init'ed:   $(Int.(isinit))"
    # @info "  infr'dims: $(infdim)"
  end
  return nothing
end

#
