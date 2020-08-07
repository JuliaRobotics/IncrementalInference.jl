
export approxConvCircular

"""
    $SIGNATURES

Construct (new) subgraph and draw the subgraph associated with clique `frontalSym::Symbol`.

Notes
- See `drawGraphCliq`/`writeGraphPdf` for details on keyword options.

Related

drawGraphCliq, spyCliqMat, drawTree, buildCliqSubgraphUp, buildSubgraphFromLabels!
"""
function drawCliqSubgraphUpMocking(fgl::G,
                                   treel::AbstractBayesTree,
                                   frontalSym::Symbol;
                                   show::Bool=true,
                                   filepath::String="/tmp/cliq_sfg.pdf",
                                   engine::AS1="sfdp",
                                   viewerapp::AS2="evince"  )::Nothing where {G <: AbstractDFG, AS1 <: AbstractString, AS2 <: AbstractString}
  #
  sfg = buildCliqSubgraphUp(fgl, treel, frontalSym)
  writeGraphPdf(sfg, show=show, viewerapp=viewerapp, engine=engine, filepath=filepath)
  nothing
end


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

Related

drawGraphCliq, drawTree, printCliqSummary, spyCliqMat
"""
function drawGraph(fgl::AbstractDFG;
                   viewerapp::AbstractString="evince",
                   filepath::AbstractString="/tmp/caesar/random/fg.pdf",
                   engine::AbstractString="neato", #sfdp
                   show::Bool=true )
  #
  mkpath(joinpath( "/", (split(filepath, '/')[1:(end-1)])...) )

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

"""
    $SIGNATURES

Draw the factor graph from a clique state machine history at a particular step as pdf and show.

Related

drawCliqSubgraphUpMocking, drawGraph, drawTree
"""
function drawGraphCliq(hists::Dict{Int, <: Tuple},
                       step::Int,
                       tree::AbstractBayesTree,
                       frontal::Symbol;
                       show::Bool=true  )
  #
  cid = getClique(tree, frontal).index
  cfg = hists[cid][step][4].cliqSubFg
  drawGraph(cfg, show=show)
end


"""
    $SIGNATURES

Print basic statistics about a clique variables and factors.

Related

printCliqHistorySummary
"""
function printCliqSummary(dfg::G,
                          cliq::TreeClique,
                          logger=ConsoleLogger() ) where G <: AbstractDFG
  #
  frtl = getCliqFrontalVarIds(cliq)
  seps = getCliqSeparatorVarIds(cliq)
  fcts = getCliqFactorIdsAll(cliq)

  isinit = map(x->isInitialized(dfg,x), [frtl;seps])
  infdim = map(x->getVariableInferredDim(dfg, x), [frtl;seps])

  with_logger(logger) do
    @info "Clique $(cliq.index) summary:"
    @info "  num frontals:    $(length(frtl))"
    @info "  num separators:  $(length(seps))"
    @info "  num factors:     $(length(fcts))"
    @info "  num initialized: $(sum(isinit)) of $(length(isinit))"
    @info ""
    @info "  frontals:  $(frtl)"
    @info "  separator: $(seps)"
    @info "  factors:   $(fcts)"
    @info "  init'ed:   $(Int.(isinit))"
    @info "  infr'dims: $(infdim)"
  end
  nothing
end



"""
    $SIGNATURES

Print basic summary of graph to `logger=ConsoleLogger()`.
"""
function printGraphSummary(dfg::G, logger=ConsoleLogger())::Nothing where {G <: AbstractDFG}
    vars = ls(dfg)
    fcts = lsf(dfg)

    prio = lsfPriors(dfg)

    isinit = map(x->isInitialized(dfg,x), vars)
    infdim = map(x->getVariableInferredDim(dfg, x), vars)
    numedges = map(v->length(ls(dfg, v)), vars)
    numfed = map(fc->length(ls(dfg, fc)), fcts)
    vardims = map(v->getDimension(getVariable(dfg, v)), vars)
    fctdims = map(v->getDimension(getFactor(dfg, v)), fcts)
    priodims = map(v->getDimension(getFactor(dfg, v)), prio)

    with_logger(logger) do
      @info "Distributed Factor Graph summary:"
      @info "  num variables:    $(length(vars))"
      @info "  num factors:      $(length(fcts)), w/ $(length(prio)) priors"
      @info "  var initialized:  $(sum(isinit))"
      @info ""
      @info "  var num edges: min. $(minimum(numedges)) | mean $(round(Statistics.mean(numedges),digits=2)) | 90% $(round(quantile(numedges,0.9),digits=2)) | max. $(maximum(numedges))"
      @info "  fct num edges: min. $(minimum(numfed)) | mean $(round(Statistics.mean(numfed),digits=2)) | 90% $(round(quantile(numfed,0.9),digits=2)) | max. $(maximum(numfed))"
      @info "  Variable dims: min. $(minimum(vardims)) | mean $(round(Statistics.mean(vardims),digits=2)) | 90% $(round(quantile(vardims,0.9),digits=2)) | max. $(maximum(vardims))"
      @info "  Factor dims:   min. $(minimum(fctdims)) | mean $(round(Statistics.mean(fctdims),digits=2)) | 90% $(round(quantile(fctdims,0.9),digits=2)) | max. $(maximum(fctdims))"
      @info "  Prior dimens:  min. $(minimum(priodims)) | mean $(round(Statistics.mean(priodims),digits=2)) | 90% $(round(quantile(priodims,0.9),digits=2)) | max. $(maximum(priodims))"
      @info "  var infr'dims: min. $(minimum(infdim)) | mean $(round(Statistics.mean(infdim),digits=2)) | 90% $(round(quantile(infdim,0.9),digits=2)) | max. $(maximum(infdim))"
    end
    nothing
end

"""
    $SIGNATURES

Print basic summary of graph to `logger=ConsoleLogger()`.
"""
function printSummary(dfg::G, logger=ConsoleLogger()) where G <: AbstractDFG
    printGraphSummary(dfg, logger)
end


"""
    $SIGNATURES

Build an approximate density `[Y|X,DX,.]=[X|Y,DX][DX|.]` as proposed by the conditional convolution.

Notes
- Assume both are on circular manifold, `manikde!(pts, (:Circular,))`
"""
function approxConvCircular(pX::BallTreeDensity, pDX::BallTreeDensity; N::Int=100)
  #

  # building basic factor graph
  tfg = initfg()
  addVariable!(tfg, :s1, Sphere1)
  addVariable!(tfg, :s2, Sphere1)
  addFactor!(tfg, [:s1;:s2], Sphere1Sphere1(pDX), graphinit=false)
  initManual!(tfg,:s1, pX)

  # solve for outgoing proposal value
  approxConv(tfg,:s1s2f1,:s2)
end

function approxConvCircular(pX::BallTreeDensity, pDX::SamplableBelief; N::Int=100)
  pts = reshape(rand(pDX, N), 1, :)
  pC = manikde!(pts, Sphere1)
  approxConvCircular(pX, pC)
end


function approxConvCircular(pX::SamplableBelief, pDX::BallTreeDensity; N::Int=100)
  pts = reshape(rand(pX, N), 1, :)
  pC = manikde!(pts, Sphere1)
  approxConvCircular(pC, pDX)
end



#
