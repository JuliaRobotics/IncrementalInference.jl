
"""
    $SIGNATURES

Construct (new) subgraph and draw the subgraph associated with clique `frontalSym::Symbol`.

Notes
- See `drawGraphCliq`/`writeGraphPdf` for details on keyword options.

Related

drawGraphCliq, spyCliqMat, drawTree, buildCliqSubgraphUp, buildSubgraphFromLabels
"""
function drawCliqSubgraphUp(fgl::FactorGraph,
                            treel::BayesTree,
                            frontalSym::Symbol;
                            show::Bool=true,
                            filepath::String="/tmp/cliq_sfg.pdf",
                            engine::AS1="sfdp",
                            viewerapp::AS2="evince"  )::FactorGraph where {AS1 <: AbstractString, AS2 <: AbstractString}
  #
  sfg = buildCliqSubgraphUp(fgl, treel, frontalSym)
  writeGraphPdf(sfg, show=show, viewerapp=viewerapp, engine=engine, filepath=filepath)
  return sfg
end


"""
    $SIGNATURES

Draw the factor graph `fgl` as pdf and show.

Related

drawGraphCliq, drawTree
"""
function drawGraph(fgl::G; show::Bool=true) where {G <: AbstractDFG}
  writeGraphPdf(fgl, show=show)
end

"""
    $SIGNATURES

Draw the factor graph from a clique state machine history at a particular step as pdf and show.

Related

drawCliqSubgraphUp, drawGraph, drawTree
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
