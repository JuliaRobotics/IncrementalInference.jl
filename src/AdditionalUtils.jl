
"""
    $SIGNATURES

Construct (new) subgraph and draw the subgraph associated with clique `frontalSym::Symbol`.

Notes
- See `drawGraphCliq`/`writeGraphPdf` for details on keyword options.

Related

drawGraphCliq, spyCliqMat, drawTree, buildCliqSubgraphUp, buildSubgraphFromLabels
"""
function drawCliqSubgraphUpMocking(fgl::G,
                                   treel::BayesTree,
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

Draw the factor graph `fgl` as pdf and show.

Related

drawGraphCliq, drawTree
"""
function drawGraph(fgl::G;
                   viewerapp::String="evince",
                   filepath::AS="/tmp/caesar/random/fg.pdf",
                   engine::AS="neato", #sfdp
                   show::Bool=true ) where {G <: AbstractDFG, AS <: AbstractString}
  writeGraphPdf(fgl, filepath=filepath, show=show, viewerapp=viewerapp, engine=engine)
end

"""
    $SIGNATURES

Draw the factor graph from a clique state machine history at a particular step as pdf and show.

Related

drawCliqSubgraphUpMocking, drawGraph, drawTree
"""
function drawGraphCliq(hists::Dict{Int, <: Tuple},
                       step::Int,
                       tree::BayesTree,
                       frontal::Symbol;
                       show::Bool=true  )
  #
  cid = getCliq(tree, frontal).index
  cfg = hists[cid][step][4].cliqSubFg
  writeGraphPdf(cfg, show=show)
end


"""
    $SIGNATURES

Print basic statistics about a clique variables and factors.

Related

printCliqHistorySummary
"""
function printCliqSummary(dfg::G,
                          cliq::Graphs.ExVertex,
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

#
