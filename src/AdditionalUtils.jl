
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
                                   filepath::String="/tmp/caesar/random/cliq_sfg.dot",
                                   engine::AS1="sfdp",
                                   viewerapp::AS2="xdot"  ) where {G <: AbstractDFG, AS1 <: AbstractString, AS2 <: AbstractString}
  #
  sfg = buildCliqSubgraphUp(fgl, treel, frontalSym)
  drawGraph(sfg, show=show, viewerapp=viewerapp, engine=engine, filepath=filepath)
  nothing
end


"""
    $SIGNATURES

Draw and show the factor graph `<:AbstractDFG` via system graphviz and xdot app.

Notes
- Requires system install on Linux of `sudo apt-get install xdot`
- Should not be calling outside programs.
- Need long term solution
- DFG's `toDotFile` a better solution -- view with `xdot` application.
- also try `engine={"sfdp","fdp","dot","twopi","circo","neato"}`

Notes:
- Calls external system application `xdot` to read the `.dot` file format
  - ```toDot(fg,file=...); @async run(`xdot file.dot`)```

Related

drawGraphCliq, [`drawTree`](@ref), printCliqSummary, spyCliqMat
"""
function drawGraph( fgl::AbstractDFG;
                    viewerapp::AbstractString="xdot",
                    filepath::AbstractString="/tmp/caesar/random/fg.dot",
                    engine::AbstractString="neato", #sfdp
                    show::Bool=true )
  #
  mkpath(dirname(filepath))
  #   mkpath(joinpath( "/", (split(filepath, '/')[1:(end-1)])...) )  

  @debug "Writing factor graph file"
  fext = split(filepath, '.')[end]
  fpwoext = filepath[1:(end-length(fext)-1)] # split(filepath, '.')[end-1]
  dotfile = fpwoext*".dot"

  # create the dot file
  DFG.toDotFile(fgl, dotfile)

  try
    # run(`$(engine) $(dotfile) -T$(fext) -o $(filepath)`)
    show ? (@async run(`$(viewerapp) $(dotfile)`)) : nothing
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
function drawGraphCliq( hists::Dict{Int, <: Tuple},
                        step::Int,
                        tree::AbstractBayesTree,
                        frontal::Symbol;
                        show::Bool=true  )
  #
  cid = getId(getClique(tree, frontal))
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
function approxConvCircular(pX::Union{<:BallTreeDensity,<:ManifoldKernelDensity}, 
                            pDX::Union{<:BallTreeDensity,<:ManifoldKernelDensity}; N::Int=100)
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

function approxConvCircular(pX::Union{<:BallTreeDensity,<:ManifoldKernelDensity}, 
                            pDX::SamplableBelief; N::Int=100)
  #
  pts = reshape(rand(pDX, N), 1, :)
  pC = manikde!(pts, Sphere1)
  approxConvCircular(pX, pC)
end


function approxConvCircular(pX::SamplableBelief, 
                            pDX::Union{<:BallTreeDensity,<:ManifoldKernelDensity}; N::Int=100)
  #
    pts = reshape(rand(pX, N), 1, :)
  pC = manikde!(pts, Sphere1)
  approxConvCircular(pC, pDX)
end



#
